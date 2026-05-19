"""
gene_translator.py

MGI-backed mouse gene symbol translator with optional human ortholog mapping.

Bridges symbol vintages between datasets aligned to different Ensembl/GENCODE
releases by routing through MGI canonical IDs. Backbone is MGI MRK_List1.rpt
(includes withdrawn symbols with forward replacement pointers).

Human ortholog mapping uses MGI HOM_MouseHumanSequence.rpt, which provides
MGI ID -> human HGNC symbol + Entrez ID via HomoloGene/HGNC.

Usage:
    import urllib.request
    urllib.request.urlretrieve(
        'https://raw.githubusercontent.com/dtabuena/Resources/main/Gene_Translator/gene_translator.py',
        'gene_translator.py'
    )
    %run gene_translator.py

    trans = GeneTranslator(cache_dir='./gene_translator')

    # Mouse-only (original behavior):
    trans.build(force=False)
    out = trans.translate('SEPT7', on_ambiguous='best')          # -> 'Septin7'
    df  = trans.translate_many(adata.var_names, on_ambiguous='best')

    # With human ortholog support:
    trans.build(force=False, include_orthologs=True)
    df  = trans.translate_many(adata.var_names, on_ambiguous='best', target='human')
    # output column contains human HGNC symbol; route includes 'no_ortholog' for misses

Cache format: HDF5 via h5py (no PyTables/pyarrow dependencies).
"""
import os
import json
import time
import warnings
import urllib.request
import urllib.error
import pandas as pd
import numpy as np
import h5py


class GeneTranslator:

    MGI_BASE_URL = 'https://www.informatics.jax.org/downloads/reports'
    MRK_LIST1_FILE = 'MRK_List1.rpt'
    HOM_FILE = 'HOM_MouseHumanSequence.rpt'

    MRK_LIST1_COLUMNS = [
        'mgi_id',
        'chromosome',
        'cm_position',
        'genome_start',
        'genome_end',
        'strand',
        'symbol',
        'status',
        'name',
        'marker_type',
        'feature_types',
        'synonyms',
        'mgi_id_current',
        'symbol_current',
    ]

    def __init__(self, cache_dir):
        self.cache_dir  = cache_dir
        os.makedirs(self.cache_dir, exist_ok=True)
        self.cache_path   = os.path.join(self.cache_dir, 'gene_translator_cache.h5')
        self.meta_path    = os.path.join(self.cache_dir, 'meta.json')
        self.raw_path     = os.path.join(self.cache_dir, self.MRK_LIST1_FILE)
        self.hom_raw_path = os.path.join(self.cache_dir, self.HOM_FILE)

        self.master   = None
        self.lookup   = None
        self.meta     = None
        self.orthologs = None   # MGI ID -> human symbol/entrez; loaded on demand

    def is_built(self):
        return os.path.exists(self.cache_path) and os.path.exists(self.meta_path)

    @staticmethod
    def _df_to_h5_group(group, df):
        """Write a DataFrame to an h5 group as one dataset per column."""
        for col in df.columns:
            vals = df[col].values
            if vals.dtype == object:
                if len(vals) and isinstance(vals[0], list):
                    encoded = np.array(['|'.join(x) if isinstance(x, list) else ''
                                        for x in vals], dtype=object)
                    group.create_dataset(col, data=encoded.astype('S'),
                                         compression='gzip')
                    group[col].attrs['was_list'] = True
                else:
                    encoded = np.array([('' if pd.isna(x) else str(x)) for x in vals],
                                       dtype=object)
                    group.create_dataset(col, data=encoded.astype('S'),
                                         compression='gzip')
                    group[col].attrs['was_list'] = False
            else:
                group.create_dataset(col, data=vals, compression='gzip')
                group[col].attrs['was_list'] = False
        group.attrs['columns'] = list(df.columns)
        group.attrs['n_rows']  = len(df)

    @staticmethod
    def _h5_group_to_df(group):
        """Read an h5 group written by _df_to_h5_group back into a DataFrame."""
        cols = list(group.attrs['columns'])
        data = {}
        for col in cols:
            arr = group[col][:]
            if arr.dtype.kind == 'S':
                arr = np.array([x.decode('utf-8') for x in arr], dtype=object)
            if group[col].attrs.get('was_list', False):
                arr = np.array([x.split('|') if x else [] for x in arr], dtype=object)
            data[col] = arr
        return pd.DataFrame(data, columns=cols)

    def build(self, force=False, include_orthologs=True):
        if self.is_built() and not force:
            print(f'[build] cache exists at {self.cache_dir}, skipping (force=False)')
            self.load()
            # If orthologs requested but not yet cached, build them now
            if include_orthologs and not self._orthologs_cached():
                self._build_orthologs()
            return self

        print(f'[build] downloading {self.MRK_LIST1_FILE}...')
        url = f'{self.MGI_BASE_URL}/{self.MRK_LIST1_FILE}'
        try:
            urllib.request.urlretrieve(url, self.raw_path)
        except urllib.error.URLError as exc:
            raise RuntimeError(f'MGI download failed from {url}: {exc!r}')

        size_mb = os.path.getsize(self.raw_path) / 1e6
        print(f'[build] downloaded {size_mb:.1f} MB')

        print('[build] parsing MRK_List1.rpt...')
        raw = pd.read_csv(
            self.raw_path,
            sep='\t',
            header=0,
            names=self.MRK_LIST1_COLUMNS,
            dtype=str,
            na_values=[''],
            keep_default_na=False,
        )
        print(f'[build] parsed {len(raw):,} rows')

        # ---------- vectorized master construction ----------
        is_official  = raw['status'] == 'O'
        is_withdrawn = raw['status'] == 'W'

        canonical_id  = raw['mgi_id'].where(is_official, raw['mgi_id_current'])
        canonical_sym = raw['symbol'].where(is_official, raw['symbol_current'])

        keep_mask = is_official | (
            is_withdrawn & canonical_id.notna() & (canonical_id != '')
        )
        kept = pd.DataFrame({
            'mgi_id':         canonical_id[keep_mask].values,
            'symbol_current': canonical_sym[keep_mask].values,
            'symbol_seen':    raw.loc[keep_mask, 'symbol'].values,
            'synonyms':       raw.loc[keep_mask, 'synonyms'].fillna('').values,
            'marker_type':    raw.loc[keep_mask, 'marker_type'].values,
        })
        print(f'[build] kept {len(kept):,} rows after status filter')

        kept['synonyms_split'] = kept['synonyms'].str.split('|')

        agg = kept.groupby('mgi_id', sort=False).agg(
            symbol_current=('symbol_current', 'first'),
            marker_type=('marker_type', 'first'),
            symbols_seen=('symbol_seen', list),
            synonyms_lists=('synonyms_split', list),
        ).reset_index()

        def _flatten(symbols_seen, synonyms_lists):
            combined = set(s for s in symbols_seen if s and not pd.isna(s))
            for sl in synonyms_lists:
                if isinstance(sl, list):
                    combined.update(s for s in sl if s and not pd.isna(s))
            return sorted(combined)

        agg['all_symbols'] = [
            _flatten(s, sl)
            for s, sl in zip(agg['symbols_seen'], agg['synonyms_lists'])
        ]
        agg = agg.drop(columns=['symbols_seen', 'synonyms_lists'])
        self.master = agg

        # ---------- vectorized lookup construction ----------
        current_edges = pd.DataFrame({
            'symbol_upper': agg['symbol_current'].str.upper(),
            'mgi_id':       agg['mgi_id'],
            'kind':         'current',
        }).dropna(subset=['symbol_upper'])

        synonym_edges = agg[['mgi_id', 'symbol_current', 'all_symbols']].explode(
            'all_symbols'
        )
        synonym_edges = synonym_edges[
            synonym_edges['all_symbols'].notna()
            & (synonym_edges['all_symbols'] != synonym_edges['symbol_current'])
        ]
        synonym_edges = pd.DataFrame({
            'symbol_upper': synonym_edges['all_symbols'].str.upper(),
            'mgi_id':       synonym_edges['mgi_id'],
            'kind':         'synonym',
        })

        lookup_df = pd.concat([current_edges, synonym_edges], ignore_index=True)
        lookup_df = lookup_df.drop_duplicates(
            subset=['symbol_upper', 'mgi_id', 'kind']
        ).reset_index(drop=True)
        self.lookup = lookup_df

        # ---------- persist via h5py ----------
        with h5py.File(self.cache_path, 'w') as f:
            self._df_to_h5_group(f.create_group('master'), self.master)
            self._df_to_h5_group(f.create_group('lookup'), self.lookup)

        self.meta = {
            'built_at':       time.strftime('%Y-%m-%d %H:%M:%S'),
            'source_url':     url,
            'raw_size_bytes': os.path.getsize(self.raw_path),
            'n_loci':         int(len(self.master)),
            'n_lookups':      int(len(self.lookup)),
        }
        with open(self.meta_path, 'w') as f:
            json.dump(self.meta, f, indent=2)

        print(f'[build] master: {len(self.master):,} loci')
        print(f'[build] lookup: {len(self.lookup):,} symbol edges')
        print(f'[build] cache written to {self.cache_path}')

        if include_orthologs:
            self._build_orthologs()

        return self

    def load(self):
        if not self.is_built():
            raise RuntimeError(
                f'translator not built; call build() first '
                f'(cache_dir={self.cache_dir})'
            )
        with h5py.File(self.cache_path, 'r') as f:
            self.master = self._h5_group_to_df(f['master'])
            self.lookup = self._h5_group_to_df(f['lookup'])
            if 'orthologs' in f:
                self.orthologs = self._h5_group_to_df(f['orthologs'])
        with open(self.meta_path, 'r') as f:
            self.meta = json.load(f)
        self._build_indices()
        has_orth = self.orthologs is not None
        print(f'[load] {len(self.master):,} loci, {len(self.lookup):,} symbol edges '
              f'(built {self.meta["built_at"]})'
              + (f', {len(self.orthologs):,} mouse-human ortholog pairs' if has_orth else ''))
        return self

    def _orthologs_cached(self):
        """Check whether the ortholog table is present in the HDF5 cache."""
        if not os.path.exists(self.cache_path):
            return False
        with h5py.File(self.cache_path, 'r') as f:
            return 'orthologs' in f

    def _build_orthologs(self):
        """
        Download HOM_MouseHumanSequence.rpt and build MGI ID -> human symbol/entrez
        mapping. Appends an 'orthologs' group to the existing HDF5 cache.

        HOM file columns (tab-delimited, no header row name — first row IS header):
            HomoloGene ID, Common Organism Name, NCBI Taxon ID, Symbol,
            EntrezGene ID, Mouse MGI ID, HGNC ID, OMIM Gene ID, Genetic Location,
            Genomic Coordinates (mouse: chr|start|end  human: chr|start|end)
        """
        url = f'{self.MGI_BASE_URL}/{self.HOM_FILE}'
        print(f'[orthologs] downloading {self.HOM_FILE}...')
        try:
            urllib.request.urlretrieve(url, self.hom_raw_path)
        except urllib.error.URLError as exc:
            raise RuntimeError(f'MGI ortholog download failed from {url}: {exc!r}')

        size_mb = os.path.getsize(self.hom_raw_path) / 1e6
        print(f'[orthologs] downloaded {size_mb:.1f} MB')

        hom = pd.read_csv(self.hom_raw_path, sep='\t', dtype=str, keep_default_na=False)

        # Normalise column names: lowercase + underscores
        hom.columns = (
            hom.columns.str.strip()
                       .str.lower()
                       .str.replace(r'[\s/\(\)]+', '_', regex=True)
                       .str.replace(r'_+$', '', regex=True)
        )

        # Identify organism column (handles slight naming variations across releases)
        org_col = next(
            (c for c in hom.columns if 'organism' in c or 'common' in c),
            None
        )
        if org_col is None:
            raise RuntimeError(
                f'Cannot identify organism column in HOM file. '
                f'Columns found: {list(hom.columns)}'
            )

        sym_col = next((c for c in hom.columns if c == 'symbol'), None)
        if sym_col is None:
            raise RuntimeError(
                f'Cannot identify symbol column in HOM file. '
                f'Columns found: {list(hom.columns)}'
            )

        entrez_col = next(
            (c for c in hom.columns if 'entrez' in c or 'entrezgene' in c),
            None
        )
        homologene_col = next(
            (c for c in hom.columns if 'homologene' in c or 'homolo' in c),
            None
        )
        if homologene_col is None:
            raise RuntimeError(
                f'Cannot identify HomoloGene ID column in HOM file. '
                f'Columns found: {list(hom.columns)}'
            )

        mouse_rows = hom[hom[org_col].str.lower().str.contains('mouse', na=False)].copy()
        human_rows = hom[hom[org_col].str.lower().str.contains('human', na=False)].copy()

        print(f'[orthologs] mouse rows: {len(mouse_rows):,}  human rows: {len(human_rows):,}')

        # Build HomoloGene group ID -> human symbol + entrez
        human_lookup = human_rows[[homologene_col, sym_col]].copy()
        human_lookup.columns = ['homologene_id', 'human_symbol']
        if entrez_col:
            human_lookup['human_entrez_id'] = human_rows[entrez_col].values
        else:
            human_lookup['human_entrez_id'] = ''
        human_lookup = human_lookup.drop_duplicates(subset=['homologene_id'])

        # MGI ID col: look for 'mgi' in column name
        mgi_col = next((c for c in hom.columns if 'mgi' in c), None)

        if mgi_col is not None and mouse_rows[mgi_col].str.startswith('MGI:').any():
            # Direct MGI ID available
            mouse_lookup = mouse_rows[[homologene_col, mgi_col]].copy()
            mouse_lookup.columns = ['homologene_id', 'mgi_id']
            mouse_lookup = mouse_lookup[mouse_lookup['mgi_id'].str.startswith('MGI:', na=False)]
        else:
            # Fall back: match by mouse symbol -> our master table
            mouse_lookup = mouse_rows[[homologene_col, sym_col]].copy()
            mouse_lookup.columns = ['homologene_id', 'mouse_symbol']
            sym_to_mgi = self.master.set_index('symbol_current')['mgi_id'].to_dict()
            mouse_lookup['mgi_id'] = mouse_lookup['mouse_symbol'].map(sym_to_mgi)
            mouse_lookup = mouse_lookup.dropna(subset=['mgi_id'])
            mouse_lookup = mouse_lookup[['homologene_id', 'mgi_id']]

        # Join mouse MGI IDs to human symbols via HomoloGene group ID
        orth = mouse_lookup.merge(human_lookup, on='homologene_id', how='inner')
        orth = orth[['mgi_id', 'human_symbol', 'human_entrez_id']].drop_duplicates(
            subset=['mgi_id']
        ).reset_index(drop=True)

        print(f'[orthologs] {len(orth):,} mouse MGI ID -> human symbol pairs built')

        self.orthologs = orth

        # Append to existing HDF5 cache
        with h5py.File(self.cache_path, 'a') as f:
            if 'orthologs' in f:
                del f['orthologs']
            self._df_to_h5_group(f.create_group('orthologs'), orth)

        print(f'[orthologs] ortholog table appended to cache')

    def _ensure_orthologs(self):
        """Raise a clear error if ortholog table is not available."""
        if self.orthologs is None:
            raise RuntimeError(
                'Ortholog table not loaded. '
                'Re-run build(include_orthologs=True) to download and cache it, '
                'then reload with load().'
            )

    def _build_indices(self):
        """
        Build hash indices for O(1) lookups instead of O(n) linear scans.
        - lookup_by_symbol: symbol_upper -> array of row indices in self.lookup
        - master_by_id:     mgi_id       -> row index in self.master
        """
        self.lookup_by_symbol = self.lookup.groupby('symbol_upper').indices
        self.master_by_id = pd.Series(
            np.arange(len(self.master)),
            index=self.master['mgi_id'].values,
        ).to_dict()

    def _ensure_loaded(self):
        if self.master is None or self.lookup is None:
            self.load()
        if not hasattr(self, 'lookup_by_symbol') or self.lookup_by_symbol is None:
            self._build_indices()

    def _hits_for(self, key):
        """O(1) lookup of all rows matching symbol_upper == key."""
        idx = self.lookup_by_symbol.get(key)
        if idx is None:
            return self.lookup.iloc[[]]
        return self.lookup.iloc[idx]

    def _current_for(self, mgi_id):
        """O(1) lookup of canonical symbol for an mgi_id."""
        i = self.master_by_id.get(mgi_id)
        if i is None:
            return None
        return self.master['symbol_current'].iat[i]

    def translate(self, symbol, on_ambiguous='flag'):
        """
        Translate one symbol to its current MGI canonical symbol.

        on_ambiguous:
            'flag'   - return list of all candidates if multiple matches
            'strict' - return None and warn if multiple matches
            'best'   - prefer 'current' over 'synonym' edges
        """
        self._ensure_loaded()
        if not isinstance(symbol, str) or not symbol.strip():
            return None

        key = symbol.strip().upper()
        hits = self._hits_for(key)
        if len(hits) == 0:
            return None

        unique_ids = hits['mgi_id'].unique().tolist()

        if len(unique_ids) == 1:
            return self._current_for(unique_ids[0])

        candidates = []
        for mgi_id in unique_ids:
            current = self._current_for(mgi_id)
            kinds = hits.loc[hits['mgi_id'] == mgi_id, 'kind'].tolist()
            candidates.append({
                'mgi_id':         mgi_id,
                'symbol_current': current,
                'kinds':          kinds,
            })

        if on_ambiguous == 'flag':
            return candidates
        if on_ambiguous == 'strict':
            warnings.warn(
                f'ambiguous translation for {symbol!r}: '
                f'{[c["symbol_current"] for c in candidates]}'
            )
            return None
        if on_ambiguous == 'best':
            current_hits = [c for c in candidates if 'current' in c['kinds']]
            chosen = current_hits[0] if current_hits else candidates[0]
            return chosen['symbol_current']
        raise ValueError(f'unknown on_ambiguous={on_ambiguous!r}')

    def translate_many(self, symbols, on_ambiguous='flag', target='mouse'):
        """
        Translate a list of symbols. Returns DataFrame columns:
            input, output, mgi_id, route, n_candidates, ambiguous

        Parameters
        ----------
        symbols      : iterable of str
        on_ambiguous : 'flag' | 'strict' | 'best'
            How to handle symbols that map to multiple MGI IDs.
        target       : 'mouse' | 'human'
            'mouse' (default) returns MGI canonical mouse symbol in output column.
            'human' maps through the ortholog table to return human HGNC symbol.
            When target='human', two extra columns are added:
                human_symbol    : human HGNC symbol (same as output, for clarity)
                human_entrez_id : human Entrez gene ID
            Genes with no human ortholog get route='no_ortholog', output=None.

        route values (mouse): identity, case_only, synonym, ambiguous, unmapped
        route values (human): above + no_ortholog
        """
        self._ensure_loaded()

        if target not in ('mouse', 'human'):
            raise ValueError(f"target must be 'mouse' or 'human', got {target!r}")

        if target == 'human':
            self._ensure_orthologs()

        symbols = list(symbols)
        results = []
        for s in symbols:
            row = {
                'input':        s,
                'output':       None,
                'mgi_id':       None,
                'route':        'unmapped',
                'n_candidates': 0,
                'ambiguous':    False,
            }
            if not isinstance(s, str) or not s.strip():
                results.append(row)
                continue

            key = s.strip().upper()
            hits = self._hits_for(key)
            unique_ids = hits['mgi_id'].unique().tolist()
            row['n_candidates'] = len(unique_ids)

            if len(unique_ids) == 0:
                results.append(row)
                continue

            if len(unique_ids) > 1:
                row['ambiguous'] = True
                if on_ambiguous == 'strict':
                    row['route'] = 'ambiguous'
                    results.append(row)
                    continue
                if on_ambiguous == 'best':
                    current_hits = hits[hits['kind'] == 'current']
                    chosen_mgi = (current_hits['mgi_id'].iloc[0]
                                  if len(current_hits) > 0
                                  else hits['mgi_id'].iloc[0])
                else:
                    row['route'] = 'ambiguous'
                    chosen_mgi = hits['mgi_id'].iloc[0]
            else:
                chosen_mgi = unique_ids[0]

            current_symbol = self._current_for(chosen_mgi)
            row['mgi_id'] = chosen_mgi
            row['output'] = current_symbol

            if row['route'] != 'ambiguous':
                if current_symbol == s:
                    row['route'] = 'identity'
                elif current_symbol is not None and current_symbol.upper() == s.upper():
                    row['route'] = 'case_only'
                else:
                    kinds = hits.loc[hits['mgi_id'] == chosen_mgi, 'kind'].tolist()
                    row['route'] = 'case_only' if 'current' in kinds else 'synonym'

            results.append(row)

        df = pd.DataFrame(results)

        if target == 'human':
            # Build MGI ID -> human symbol/entrez lookup dict from ortholog table
            orth_symbol = self.orthologs.set_index('mgi_id')['human_symbol'].to_dict()
            orth_entrez  = self.orthologs.set_index('mgi_id')['human_entrez_id'].to_dict()

            human_symbols = df['mgi_id'].map(orth_symbol)
            human_entrez  = df['mgi_id'].map(orth_entrez)

            # Mark genes that resolved to mouse but have no human ortholog
            mouse_resolved = df['route'].isin(['identity', 'case_only', 'synonym', 'ambiguous'])
            no_ortholog    = mouse_resolved & human_symbols.isna()
            df.loc[no_ortholog, 'route'] = 'no_ortholog'

            # output column = human symbol (None if unmapped or no ortholog)
            df['output']          = human_symbols.where(~no_ortholog & mouse_resolved, other=None)
            df['human_symbol']    = df['output']
            df['human_entrez_id'] = human_entrez.where(~no_ortholog & mouse_resolved, other=None)

        return df

    # ---------- ANNDATA HELPERS ----------

    def canonicalize_anndata(self, adata, uppercase=True, on_ambiguous='best',
                             require_counts=True, verbose=True):
        """
        Replace adata.var_names with canonical MGI symbols (optionally uppercased).
        Sum-merges any rows whose canonical names collide. Returns a NEW AnnData
        (the input is not mutated).

        IMPORTANT: sum-merging is only valid on raw counts. By default this
        function refuses to operate on data that doesn't look like integer counts
        (set require_counts=False to override, but understand the math first).

        Parameters
        ----------
        adata           AnnData
            Must contain raw counts in adata.X.
        uppercase       bool
            If True (default), target names are UPPERCASE canonical (e.g. APOE).
            If False, target names are MGI Title-case (e.g. Apoe).
        on_ambiguous    'best' | 'flag' | 'strict'
            Passed through to translate_many for multi-MGI-ID symbols.
        require_counts  bool
            If True (default), abort if X doesn't look like integer counts.
        verbose         bool
            Print build progress.

        Returns
        -------
        new_adata : AnnData
            Same n_obs, possibly fewer n_vars (after collision merges).
            new_adata.var has columns:
                symbol_original  pipe-joined original symbol(s) per merged row
                mgi_id           canonical MGI ID (or None for unmapped)
                route            pipe-joined translation route(s)
                n_merged         number of source rows merged into this row
                unmapped         True if row's name is original (not translated)
            Original adata.obs is preserved as-is.
            adata.layers, .obsm, .uns are NOT carried over (they may be invalid
            after row aggregation; recompute downstream).
        """
        import scipy.sparse as _sp
        import anndata as _ad

        self._ensure_loaded()

        # Validate count-like input
        X = adata.X
        if require_counts:
            sample = X[:min(1000, X.shape[0]), :].toarray() if _sp.issparse(X) else X[:min(1000, X.shape[0]), :]
            if not np.allclose(sample, np.round(sample)):
                raise ValueError(
                    'adata.X does not look like integer counts. '
                    'Sum-merge is only mathematically valid on raw counts. '
                    'Pass require_counts=False to override (use at your own risk).'
                )

        # Translate
        if verbose:
            print(f'[canonicalize] translating {adata.n_vars:,} symbols...')
        tr = self.translate_many(adata.var_names, on_ambiguous=on_ambiguous)
        canonical = tr['output'].values
        unmapped_mask = pd.isna(canonical)

        # Build target names
        target = np.where(unmapped_mask, adata.var_names.values, canonical)
        if uppercase:
            target = np.array([str(t).upper() for t in target])
        else:
            target = np.array([str(t) for t in target])

        # Determine collisions
        unique_targets, inverse = np.unique(target, return_inverse=True)
        n_in = len(target)
        n_out = len(unique_targets)
        n_collisions = n_in - n_out

        if verbose:
            print(f'[canonicalize] translated: {(~unmapped_mask).sum():,}, '
                  f'unmapped: {unmapped_mask.sum():,}')
            print(f'[canonicalize] {n_in:,} input -> {n_out:,} output '
                  f'({n_collisions} rows merged)')

        # Sparse aggregation: M @ agg.T sums input cols within each target group
        if not _sp.issparse(X):
            X = _sp.csr_matrix(X)

        agg = _sp.csr_matrix(
            (np.ones(n_in, dtype=np.float32),
             (inverse, np.arange(n_in))),
            shape=(n_out, n_in),
        )
        X_out = (X @ agg.T).tocsr().astype(np.float32)

        # Sanity: total counts conserved (within float32 precision)
        total_in = float(X.sum())
        total_out = float(X_out.sum())
        rel_err = abs(total_in - total_out) / (total_in + 1e-12)
        if rel_err > 1e-5:
            raise RuntimeError(
                f'count conservation broken: {total_in:,.0f} -> {total_out:,.0f} '
                f'(rel err {rel_err:.2e})'
            )

        # Build var dataframe
        new_var_rows = []
        var_names_arr = adata.var_names.values
        mgi_ids_arr = tr['mgi_id'].values
        routes_arr = tr['route'].values
        for j in range(n_out):
            in_idx = np.where(inverse == j)[0]
            originals = var_names_arr[in_idx].tolist()
            mgi_ids = mgi_ids_arr[in_idx].tolist()
            routes = routes_arr[in_idx].tolist()
            new_var_rows.append({
                'symbol_original': '|'.join(str(o) for o in originals),
                'mgi_id':          mgi_ids[0] if mgi_ids[0] is not None else None,
                'route':           '|'.join(str(r) for r in routes),
                'n_merged':        len(in_idx),
                'unmapped':        bool(unmapped_mask[in_idx[0]]),
            })

        new_var = pd.DataFrame(
            new_var_rows,
            index=pd.Index(unique_targets, name='symbol'),
        )

        new_adata = _ad.AnnData(
            X=X_out,
            obs=adata.obs.copy(),
            var=new_var,
        )

        if verbose:
            n_multi = int((new_var['n_merged'] > 1).sum())
            n_unmapped = int(new_var['unmapped'].sum())
            print(f'[canonicalize] result: {new_adata.shape}  '
                  f'multi-source rows: {n_multi}  unmapped: {n_unmapped}')

        return new_adata
