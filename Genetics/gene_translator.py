"""
gene_translator.py

MGI-backed mouse gene symbol translator.

Bridges symbol vintages between datasets aligned to different Ensembl/GENCODE
releases by routing through MGI canonical IDs. Backbone is MGI MRK_List1.rpt
(includes withdrawn symbols with forward replacement pointers).

Usage:
    import urllib.request
    urllib.request.urlretrieve(
        'https://raw.githubusercontent.com/dtabuena/Resources/main/Gene_Translator/gene_translator.py',
        'gene_translator.py'
    )
    %run gene_translator.py

    trans = GeneTranslator(cache_dir='./gene_translator')
    trans.build(force=False)        # downloads ~86 MB once; instant on subsequent calls
    out = trans.translate('SEPT7', on_ambiguous='best')          # -> 'Septin7'
    df  = trans.translate_many(adata.var_names, on_ambiguous='best')

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
        self.cache_dir = cache_dir
        os.makedirs(self.cache_dir, exist_ok=True)
        self.cache_path = os.path.join(self.cache_dir, 'gene_translator_cache.h5')
        self.meta_path  = os.path.join(self.cache_dir, 'meta.json')
        self.raw_path   = os.path.join(self.cache_dir, self.MRK_LIST1_FILE)

        self.master = None
        self.lookup = None
        self.meta   = None

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

    def build(self, force=False):
        if self.is_built() and not force:
            print(f'[build] cache exists at {self.cache_dir}, skipping (force=False)')
            return self.load()

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
        with open(self.meta_path, 'r') as f:
            self.meta = json.load(f)
        self._build_indices()
        print(f'[load] {len(self.master):,} loci, {len(self.lookup):,} symbol edges '
              f'(built {self.meta["built_at"]})')
        return self

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

    def translate_many(self, symbols, on_ambiguous='flag'):
        """
        Translate a list of symbols. Returns DataFrame columns:
            input, output, mgi_id, route, n_candidates, ambiguous

        route values: identity, case_only, synonym, ambiguous, unmapped
        """
        self._ensure_loaded()

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

        return pd.DataFrame(results)
