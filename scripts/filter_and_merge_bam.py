#!/usr/bin/env python3
"""Filter per-sample BAM files using transcript assignments and merge retained reads."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Set

import pysam

from common import ensure_dir, load_samples_tsv, open_maybe_gzip, resolve_sample_dirs


def load_valid_read_ids(tsv_path: str) -> Set[str]:
    valid: Set[str] = set()
    with open_maybe_gzip(tsv_path, 'rt') as handle:
        reader = csv.reader(handle, delimiter='\t')
        for row in reader:
            if not row or row[0].startswith('#'):
                continue
            if len(row) < 2 or row[1] == '*':
                continue
            valid.add(row[0])
    return valid


def process_sample(sample_id: str, bam_path: str, transcript_tsv: str, out_dir: str) -> Dict[str, str | int]:
    valid_ids = load_valid_read_ids(transcript_tsv)
    out_bam = str(Path(out_dir) / f'{sample_id}_filtered.bam')
    reads_in = 0
    reads_out = 0
    with pysam.AlignmentFile(bam_path, 'rb') as in_bam, pysam.AlignmentFile(out_bam, 'wb', template=in_bam) as out_handle:
        for read in in_bam.fetch(until_eof=True):
            reads_in += 1
            if read.query_name in valid_ids:
                out_handle.write(read)
                reads_out += 1
    pysam.index(out_bam)
    return {
        'sample_id': sample_id,
        'reads_in_bam': reads_in,
        'reads_written': reads_out,
        'valid_read_ids': len(valid_ids),
        'output_bam': out_bam,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description='Filter and merge BAM files.')
    parser.add_argument('--samples', required=True)
    parser.add_argument('--project-dir', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--bam-name', default='{sample_id}_sorted.bam')
    parser.add_argument('--transcript-reads-name', default='OUT.transcript_model_reads.tsv.gz')
    parser.add_argument('--merged-name', default='all_samples_merged.bam')
    args = parser.parse_args()

    samples = resolve_sample_dirs(args.project_dir, load_samples_tsv(args.samples))
    output_dir = ensure_dir(args.output_dir)
    results: List[Dict[str, str | int]] = []

    with ProcessPoolExecutor(max_workers=args.threads) as pool:
        futures = {}
        for sample in samples:
            sample_id = sample['sample_id']
            bam_name = args.bam_name.format(sample_id=sample_id)
            futures[pool.submit(
                process_sample,
                sample_id,
                str(Path(sample['sample_dir']) / bam_name),
                str(Path(sample['out_dir']) / args.transcript_reads_name),
                str(output_dir),
            )] = sample_id
        for future in as_completed(futures):
            results.append(future.result())

    stats_path = output_dir / 'filtering_statistics.tsv'
    with open(stats_path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['sample_id', 'valid_read_ids', 'reads_in_bam', 'reads_written', 'output_bam'], delimiter='\t')
        writer.writeheader()
        writer.writerows(sorted(results, key=lambda x: x['sample_id']))

    bam_files = [str(r['output_bam']) for r in sorted(results, key=lambda x: x['sample_id'])]
    merged_bam = str(output_dir / args.merged_name)
    pysam.merge('-f', merged_bam, *bam_files)
    pysam.index(merged_bam)


if __name__ == '__main__':
    main()
