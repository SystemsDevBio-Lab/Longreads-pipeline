#!/usr/bin/env python3
"""Generate per-sample and merged gene-by-cell matrices from IsoQuant read assignments."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

from common import (
    ensure_dir,
    load_samples_tsv,
    open_maybe_gzip,
    parse_read_id,
    reindex_matrix,
    resolve_sample_dirs,
    sorted_union,
    write_lines,
    write_matrix_market_bundle,
)


def read_gene_assignments(path: str) -> List[Tuple[str, str]]:
    rows: List[Tuple[str, str]] = []
    with open_maybe_gzip(path, 'rt') as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split()
            if len(parts) < 5:
                continue
            read_id = parts[0]
            gene_id = parts[4]
            if gene_id in {'*', '.', ''}:
                continue
            rows.append((read_id, gene_id))
    return rows


def build_gene_matrix(assignments: List[Tuple[str, str]], sample_id: str, append_sample: bool = True, barcode_sep: str = '_') -> pd.DataFrame:
    bu_to_genes: Dict[str, set] = defaultdict(set)
    bu_to_barcode: Dict[str, str] = {}
    for read_id, gene_id in assignments:
        try:
            barcode, _, bu = parse_read_id(read_id)
        except ValueError:
            continue
        bu_to_genes[bu].add(gene_id)
        bu_to_barcode[bu] = barcode

    retained = []
    for bu, genes in bu_to_genes.items():
        if len(genes) == 1:
            gene_id = next(iter(genes))
            barcode = bu_to_barcode[bu]
            barcode_name = f"{sample_id}{barcode_sep}{barcode}" if append_sample else barcode
            retained.append((gene_id, barcode_name))

    if not retained:
        return pd.DataFrame()

    df = pd.DataFrame(retained, columns=['gene_id', 'cell'])
    matrix = df.groupby(['gene_id', 'cell']).size().unstack(fill_value=0)
    matrix = matrix.sort_index(axis=0).sort_index(axis=1)
    return matrix.astype('int32')


def merge_gene_matrices(matrices: List[pd.DataFrame]) -> pd.DataFrame:
    if not matrices:
        return pd.DataFrame()
    merged = pd.concat(matrices, axis=1).fillna(0).astype('int32')
    merged = merged.loc[~merged.index.duplicated(keep='first')]
    merged = merged.sort_index(axis=0).sort_index(axis=1)
    return merged


def discover_transcript_barcodes(samples: List[dict[str, str]], bu_name: str, barcode_sep: str) -> List[str]:
    columns: List[str] = []
    for sample in samples:
        sample_id = sample['sample_id']
        bu_path = Path(sample['out_dir']) / bu_name.format(sample_id=sample_id)
        if not bu_path.exists():
            continue
        df = pd.read_csv(bu_path, dtype=str)
        if 'status' in df.columns:
            df = df[df['status'] == 'unique'].copy()
        elif 'read_id' in df.columns:
            df = df[df['read_id'] != 'ambiguous'].copy()
        if 'barcode' not in df.columns:
            continue
        sample_columns = [
            f"{sample_id}{barcode_sep}{barcode}"
            for barcode in df['barcode'].dropna().astype(str).unique().tolist()
        ]
        columns.extend(sample_columns)
    return sorted_union(columns)


def main() -> None:
    parser = argparse.ArgumentParser(description='Build per-sample and merged gene matrices.')
    parser.add_argument('--samples', required=True, help='Path to samples.tsv')
    parser.add_argument('--project-dir', required=True, help='Project root directory')
    parser.add_argument('--assignment-name', default='OUT.read_assignments.tsv.gz')
    parser.add_argument('--output-name', default='gene_matrix.csv.gz')
    parser.add_argument('--merged-output', default='multisample/matrices/gene_matrix.csv.gz')
    parser.add_argument('--merged-output-uncompressed', default='multisample/matrices/gene_matrix.csv')
    parser.add_argument('--merged-mtx-prefix', default='multisample/matrices/gene_matrix')
    parser.add_argument('--barcode-sep', default='_')
    parser.add_argument('--bu-name', default='{sample_id}_bu.csv')
    parser.add_argument('--shared-barcodes-out', default='multisample/matrices/shared_barcodes.txt')
    args = parser.parse_args()

    samples = resolve_sample_dirs(args.project_dir, load_samples_tsv(args.samples))
    matrices: List[pd.DataFrame] = []
    observed_gene_columns: List[str] = []

    for sample in samples:
        sample_id = sample['sample_id']
        assignment_path = str(Path(sample['out_dir']) / args.assignment_name)
        matrix = build_gene_matrix(
            read_gene_assignments(assignment_path),
            sample_id,
            append_sample=True,
            barcode_sep=args.barcode_sep,
        )
        out_path = Path(sample['out_dir']) / args.output_name
        ensure_dir(out_path.parent)
        matrix.to_csv(out_path, compression='gzip')
        matrices.append(matrix)
        observed_gene_columns.extend(matrix.columns.tolist())

    merged = merge_gene_matrices(matrices)
    transcript_columns = discover_transcript_barcodes(samples, args.bu_name, args.barcode_sep)
    shared_columns = sorted_union(observed_gene_columns + transcript_columns)
    merged = reindex_matrix(merged, row_names=merged.index.tolist(), col_names=shared_columns)

    merged_path = Path(args.project_dir) / args.merged_output
    merged_csv_path = Path(args.project_dir) / args.merged_output_uncompressed
    merged_mtx_prefix = Path(args.project_dir) / args.merged_mtx_prefix
    shared_barcodes_path = Path(args.project_dir) / args.shared_barcodes_out

    ensure_dir(merged_path.parent)
    merged.to_csv(merged_path, compression='gzip')
    merged.to_csv(merged_csv_path)
    write_matrix_market_bundle(
        merged.to_numpy(dtype='int32', copy=False),
        merged.index.tolist(),
        merged.columns.tolist(),
        merged_mtx_prefix,
    )
    write_lines(shared_barcodes_path, merged.columns.tolist())


if __name__ == '__main__':
    main()
