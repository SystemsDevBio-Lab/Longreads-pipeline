#!/usr/bin/env python3
"""Build the final transcript-by-cell matrix from second-pass IsoQuant output and per-sample BU tables."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

import numpy as np
import pandas as pd
from scipy import sparse

from common import (
    ensure_dir,
    load_samples_tsv,
    open_maybe_gzip,
    reindex_matrix,
    resolve_sample_dirs,
    sorted_union,
    write_lines,
    write_matrix_market_bundle,
)


def load_bu_tables(samples_tsv: str, project_dir: str, barcode_sep: str, bu_name: str) -> Tuple[Dict[str, int], List[str]]:
    samples = resolve_sample_dirs(project_dir, load_samples_tsv(samples_tsv))
    all_columns: List[str] = []
    read_to_colidx: Dict[str, int] = {}
    offset = 0
    for sample in samples:
        sample_id = sample['sample_id']
        bu_path = Path(sample['out_dir']) / bu_name.format(sample_id=sample_id)
        df = pd.read_csv(bu_path, dtype=str)
        df = df[df['status'] == 'unique'].copy() if 'status' in df.columns else df[df['read_id'] != 'ambiguous'].copy()
        barcodes = df['barcode'].dropna().astype(str).unique().tolist()
        columns = [f'{sample_id}{barcode_sep}{barcode}' for barcode in barcodes]
        col_idx_map = {name: i + offset for i, name in enumerate(columns)}
        barcode_to_global = {barcode: col_idx_map[f'{sample_id}{barcode_sep}{barcode}'] for barcode in barcodes}
        for _, row in df.iterrows():
            read_id = row['read_id']
            if read_id == 'ambiguous':
                continue
            read_to_colidx[read_id] = barcode_to_global[str(row['barcode'])]
        all_columns.extend(columns)
        offset += len(columns)
    return read_to_colidx, all_columns


def find_unique_reads(map_tsv: str) -> Set[str]:
    seen_once: Set[str] = set()
    seen_multi: Set[str] = set()
    with open_maybe_gzip(map_tsv, 'rt') as handle:
        reader = csv.reader(handle, delimiter='\t')
        for row in reader:
            if not row or row[0].startswith('#'):
                continue
            read_id = row[0]
            if read_id in seen_multi:
                continue
            if read_id in seen_once:
                seen_once.remove(read_id)
                seen_multi.add(read_id)
            else:
                seen_once.add(read_id)
    return seen_once


def build_sparse_matrix(map_tsv: str, read_to_colidx: Dict[str, int], unique_reads: Set[str]) -> Tuple[sparse.csr_matrix, List[str]]:
    transcript_to_rowidx: Dict[str, int] = {}
    counts: Dict[Tuple[int, int], int] = defaultdict(int)
    with open_maybe_gzip(map_tsv, 'rt') as handle:
        reader = csv.reader(handle, delimiter='\t')
        for row in reader:
            if not row or row[0].startswith('#'):
                continue
            if len(row) < 2:
                continue
            read_id, transcript_id = row[0], row[1]
            if transcript_id == '*' or read_id not in unique_reads:
                continue
            col_idx = read_to_colidx.get(read_id)
            if col_idx is None:
                continue
            row_idx = transcript_to_rowidx.setdefault(transcript_id, len(transcript_to_rowidx))
            counts[(row_idx, col_idx)] += 1

    n_rows = len(transcript_to_rowidx)
    n_cols = max(read_to_colidx.values()) + 1 if read_to_colidx else 0
    rows, cols, data = [], [], []
    for (r, c), v in counts.items():
        rows.append(r)
        cols.append(c)
        data.append(v)
    matrix = sparse.coo_matrix(
        (np.asarray(data, dtype=np.int32), (np.asarray(rows, dtype=np.int32), np.asarray(cols, dtype=np.int32))),
        shape=(n_rows, n_cols),
        dtype=np.int32,
    ).tocsr()
    transcripts = [None] * n_rows
    for tid, idx in transcript_to_rowidx.items():
        transcripts[idx] = tid
    return matrix, transcripts


def load_shared_barcodes(shared_barcodes_path: str | None, fallback_columns: List[str]) -> List[str]:
    if not shared_barcodes_path:
        return sorted_union(fallback_columns)
    path = Path(shared_barcodes_path)
    if not path.exists():
        return sorted_union(fallback_columns)
    with open(path, 'r') as handle:
        values = [line.strip() for line in handle if line.strip()]
    return values if values else sorted_union(fallback_columns)


def main() -> None:
    parser = argparse.ArgumentParser(description='Build transcript-by-cell matrix.')
    parser.add_argument('--map-tsv', required=True)
    parser.add_argument('--samples', required=True)
    parser.add_argument('--project-dir', required=True)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--barcode-sep', default='_')
    parser.add_argument('--bu-name', default='{sample_id}_bu.csv')
    parser.add_argument('--shared-barcodes', default=None)
    parser.add_argument('--emit-nonzero-tsv', action='store_true')
    args = parser.parse_args()

    outdir = ensure_dir(args.outdir)
    read_to_colidx, columns = load_bu_tables(args.samples, args.project_dir, args.barcode_sep, args.bu_name)
    unique_reads = find_unique_reads(args.map_tsv)
    matrix, transcripts = build_sparse_matrix(args.map_tsv, read_to_colidx, unique_reads)

    shared_columns = load_shared_barcodes(args.shared_barcodes, columns)
    matrix_df = pd.DataFrame.sparse.from_spmatrix(matrix, index=transcripts, columns=columns)
    matrix_df = matrix_df.sparse.to_dense()
    matrix_df = reindex_matrix(matrix_df, row_names=transcripts, col_names=shared_columns)
    dense_values = matrix_df.to_numpy(dtype=np.int32, copy=False)
    matrix = sparse.csr_matrix(dense_values)
    transcripts = matrix_df.index.tolist()
    columns = matrix_df.columns.tolist()

    sparse.save_npz(outdir / 'matrix.npz', matrix)
    write_lines(outdir / 'transcripts.txt', transcripts)
    write_lines(outdir / 'barcodes.txt', columns)
    matrix_df.to_csv(outdir / 'transcript_matrix.csv.gz', compression='gzip')
    matrix_df.to_csv(outdir / 'transcript_matrix.csv')
    write_matrix_market_bundle(matrix, transcripts, columns, outdir / 'transcript_matrix')

    if args.emit_nonzero_tsv:
        with open_maybe_gzip(outdir / 'nonzero.tsv.gz', 'wt') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['transcript_id', 'sampleid_barcode', 'count'])
            coo = matrix.tocoo()
            for r, c, v in zip(coo.row, coo.col, coo.data):
                writer.writerow([transcripts[r], columns[c], int(v)])


if __name__ == '__main__':
    main()
