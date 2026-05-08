#!/usr/bin/env python3
"""Build representative barcode-UMI tables for transcript quantification."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List

import pandas as pd
import pysam

from common import ensure_dir, open_maybe_gzip, parse_read_id


def read_transcript_map(tsv_path: str) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with open_maybe_gzip(tsv_path, 'rt') as handle:
        reader = csv.reader(handle, delimiter='\t')
        for row in reader:
            if not row or row[0].startswith('#'):
                continue
            if len(row) < 2:
                continue
            read_id, transcript_id = row[0], row[1]
            mapping[read_id] = transcript_id
    return mapping


def add_read_lengths(read_to_transcript: Dict[str, str], bam_path: str) -> pd.DataFrame:
    records: List[dict] = []
    seen = set(read_to_transcript)
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for read in bam.fetch(until_eof=True):
            read_id = read.query_name
            if read_id not in seen:
                continue
            records.append({
                'read_id': read_id,
                'transcript_id': read_to_transcript[read_id],
                'read_length': int(read.query_length or 0),
            })
    return pd.DataFrame(records)


def build_bu_table(df: pd.DataFrame) -> pd.DataFrame:
    df = df[df['transcript_id'] != '*'].copy()
    if df.empty:
        return pd.DataFrame(columns=['bu', 'barcode', 'umi', 'read_id', 'transcript_id', 'read_length', 'status'])

    parsed = df['read_id'].apply(lambda x: parse_read_id(x))
    df[['barcode', 'umi', 'bu']] = pd.DataFrame(parsed.tolist(), index=df.index)

    reps: List[dict] = []
    for bu, group in df.groupby('bu', sort=False):
        transcripts = set(group['transcript_id'])
        barcode = group['barcode'].iloc[0]
        umi = group['umi'].iloc[0]
        if len(transcripts) == 1:
            rep = group.sort_values(['read_length', 'read_id'], ascending=[False, True]).iloc[0]
            reps.append({
                'bu': bu,
                'barcode': barcode,
                'umi': umi,
                'read_id': rep['read_id'],
                'transcript_id': rep['transcript_id'],
                'read_length': int(rep['read_length']),
                'status': 'unique',
            })
        else:
            reps.append({
                'bu': bu,
                'barcode': barcode,
                'umi': umi,
                'read_id': 'ambiguous',
                'transcript_id': 'ambiguous',
                'read_length': 'ambiguous',
                'status': 'ambiguous',
            })
    return pd.DataFrame(reps)


def main() -> None:
    parser = argparse.ArgumentParser(description='Generate representative BU tables.')
    parser.add_argument('--transcript-reads', required=True, help='OUT.transcript_model_reads.tsv.gz')
    parser.add_argument('--bam', required=True, help='Sample sorted BAM file')
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--sample-id', required=True)
    args = parser.parse_args()

    output_dir = ensure_dir(args.output_dir)
    mapping = read_transcript_map(args.transcript_reads)
    df = add_read_lengths(mapping, args.bam)
    with_length_path = output_dir / 'OUT.transcript_model_reads_with_length.tsv.gz'
    df.to_csv(with_length_path, sep='\t', index=False, compression='gzip')

    bu_df = build_bu_table(df)
    bu_path = output_dir / f'{args.sample_id}_bu.csv'
    bu_df.to_csv(bu_path, index=False)


if __name__ == '__main__':
    main()
