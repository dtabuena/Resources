"""
gene_translator.py

MGI-backed mouse gene symbol translator.

Usage:
    import urllib.request
    urllib.request.urlretrieve(
        'https://raw.githubusercontent.com/dtabuena/Resources/main/Gene_Translator/gene_translator.py',
        'gene_translator.py'
    )
    %run gene_translator.py

    trans = GeneTranslator(cache_dir='./gene_translator')
    trans.build(force=False)   # one-time, ~30 sec download

    out = trans.translate('SEPT7', on_ambiguous='best')
    df  = trans.translate_many(adata.var_names, on_ambiguous='best')

Backbone: MGI MRK_List1.rpt (https://www.informatics.jax.org/downloads/reports).
"""
import os
import json
import time
import warnings
import urllib.request
import urllib.error
import pandas as pd
import numpy as np


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
        self.master_path = os.path.join(self.cache_dir, 'mgi_master.parquet')
        self.lookup_path = os.path.join(self.cache_dir, 'symbol_lookup.parquet')
        self.meta_path   = os.path.join(self.cache_dir, 'meta.json')
        self.raw_path    = os.path.join(self.cache_dir, self.MRK_LIST1_FILE)

        self.master = None
        self.lookup = None
        self.meta   = None

    def is_built(self):
        return (os.path.exists(self.master_path) and
                os.path.exists(self.lookup_path) and
                os.path.exists(self.meta_path))

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

        master_rows = []
        for _, row in raw.iterrows():
            status = row['status']
            sym = row['symbol']
            syns = row['synonyms']
            syn_list = [] if (pd.isna(syns) or syns == '') else syns.split('|')

            if status == 'O':
                canonical_id = row['mgi_id']
                canonical_symbol = sym
            elif status == 'W':
                canonical_id = row['mgi_id_current']
                canonical_symbol = row['symbol_current']
                if pd.isna(canonical_id) or canonical_id == '':
                    continue
            else:
                continue

            master_rows.append({
                'mgi_id':         canonical_id,
                'symbol_current': canonical_symbol,
                'symbol_seen':    sym,
                'status_seen':    status,
                'synonyms':       syn_list,
                'marker_type':    row['marker_type'],
            })

        master_df = pd.DataFrame(master_rows)

        agg = master_df.groupby('mgi_id', sort=False).agg(
            symbol_current=('symbol_current', 'first'),
            marker_type=('marker_type', 'first'),
            symbols_seen=('symbol_seen', lambda x: list(x)),
            synonyms_lists=('synonyms', lambda x: list(x)),
        )

        all_synonyms = []
        for symbols_seen, synonyms_lists in zip(
            agg['symbols_seen'], agg['synonyms_lists']
        ):
            combined = set()
            combined.update(symbols_seen)
            for sl in synonyms_lists:
                combined.update(sl)
            all_synonyms.append(sorted(s for s in combined if s and not pd.isna(s)))

        agg['all_symbols'] = all_synonyms
        agg = agg.drop(columns=['symbols_seen', 'synonyms_lists']).reset_index()
        self.master = agg

        lookup_rows = []
        for _, row in agg.iterrows():
            mgi_id = row['mgi_id']
            current = row['symbol_current']
            if current and not pd.isna(current):
                lookup_rows.append({
                    'symbol_upper': current.upper(),
                    'mgi_id':       mgi_id,
                    'kind':         'current',
                })
            for s in row['all_symbols']:
                if s == current:
                    continue
                lookup_rows.append({
                    'symbol_upper': s.upper(),
                    'mgi_id':       mgi_id,
                    'kind':         'synonym',
                })

        lookup_df = pd.DataFrame(lookup_rows).drop_duplicates(
            subset=['symbol_upper', 'mgi_id', 'kind']
        ).reset_index(drop=True)
        self.lookup = lookup_df

        self.master.to_parquet(self.master_path, index=False)
        self.lookup.to_parquet(self.lookup_path, index=False)

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
        print(f'[build] cache written to {self.cache_dir}')
        return self

    def load(self):
        if not self.is_built():
            raise RuntimeError(
                f'translator not built; call build() first '
                f'(cache_dir={self.cache_dir})'
            )
        self.master = pd.read_parquet(self.master_path)
        self.lookup = pd.read_parquet(self.lookup_path)
        with open(self.meta_path, 'r') as f:
            self.meta = json.load(f)
        print(f'[load] {len(self.master):,} loci, {len(self.lookup):,} symbol edges '
              f'(built {self.meta["built_at"]})')
        return self

    def _ensure_loaded(self):
        if self.master is None or self.lookup is None:
            self.load()

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
        hits = self.lookup[self.lookup['symbol_upper'] == key]
        if len(hits) == 0:
            return None

        unique_ids = hits['mgi_id'].unique().tolist()

        if len(unique_ids) == 1:
            mgi_id = unique_ids[0]
            current = self.master.loc[
                self.master['mgi_id'] == mgi_id, 'symbol_current'
            ].iloc[0]
            return current

        candidates = []
        for mgi_id in unique_ids:
            current = self.master.loc[
                self.master['mgi_id'] == mgi_id, 'symbol_current'
            ].iloc[0]
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
            hits = self.lookup[self.lookup['symbol_upper'] == key]
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

            mrow = self.master.loc[self.master['mgi_id'] == chosen_mgi].iloc[0]
            current_symbol = mrow['symbol_current']
            row['mgi_id'] = chosen_mgi
            row['output'] = current_symbol

            if row['route'] != 'ambiguous':
                if current_symbol == s:
                    row['route'] = 'identity'
                elif current_symbol.upper() == s.upper():
                    row['route'] = 'case_only'
                else:
                    kinds = hits.loc[hits['mgi_id'] == chosen_mgi, 'kind'].tolist()
                    row['route'] = 'case_only' if 'current' in kinds else 'synonym'

            results.append(row)

        return pd.DataFrame(results)
