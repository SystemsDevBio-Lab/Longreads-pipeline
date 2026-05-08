#!/usr/bin/env python3
"""Collapse transcript GTFs across samples into a consensus annotation."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from common import ensure_dir, load_samples_tsv, resolve_sample_dirs

MIN_READS_THRESHOLD = 5
MIN_SAMPLES_FOR_NOVEL = 2


@dataclass
class Exon:
    start: int
    end: int


@dataclass
class Transcript:
    chrom: str
    strand: str
    exons: List[Exon]
    gene_id: str
    transcript_id: str
    source: str
    sample_id: str
    is_reference: bool

    @property
    def exon_count(self) -> int:
        return len(self.exons)

    @property
    def start(self) -> int:
        return self.exons[0].start

    @property
    def end(self) -> int:
        return self.exons[-1].end

    def ends(self) -> Tuple[int, int]:
        if self.strand == '+':
            return self.start, self.end
        return self.end, self.start

    def junctions(self) -> List[Tuple[int, int]]:
        return [(self.exons[i].end, self.exons[i + 1].start) for i in range(self.exon_count - 1)]


def parse_attrs(attr_str: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for field in attr_str.strip().strip(';').split(';'):
        field = field.strip()
        if not field or ' ' not in field:
            continue
        key, value = field.split(' ', 1)
        attrs[key] = value.strip().strip('"')
    return attrs


def is_reference(transcript_id: str, gene_id: str, source: str, attrs: Dict[str, str]) -> bool:
    return gene_id.startswith('ENSG') or transcript_id.startswith('ENS') or source in {'ENSEMBL', 'HAVANA'} or 'havana_transcript' in attrs


def normalize_exons(exons: List[Exon]) -> List[Exon]:
    exons = sorted(exons, key=lambda x: (x.start, x.end))
    return exons


def read_transcript_counts(matrix_csv: Optional[str], min_reads: int) -> set[str]:
    if not matrix_csv or not Path(matrix_csv).exists():
        return set()
    keep: set[str] = set()
    with open(matrix_csv, 'r', newline='') as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        for row in reader:
            if not row:
                continue
            transcript_id = row[0]
            try:
                total = sum(float(x) for x in row[1:] if x != '')
            except ValueError:
                total = 0
            if total >= min_reads:
                keep.add(transcript_id)
    return keep


def read_gtf(gtf_path: str, sample_id: str, min_reads: int, matrix_csv: Optional[str]) -> List[Transcript]:
    passing = read_transcript_counts(matrix_csv, min_reads)
    transcript_meta: Dict[str, Tuple[str, str, str, str, bool]] = {}
    exon_map: Dict[str, List[Exon]] = defaultdict(list)
    with open(gtf_path, 'r') as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs_str = line.rstrip('\n').split('\t')
            attrs = parse_attrs(attrs_str)
            tid = attrs.get('transcript_id')
            gid = attrs.get('gene_id')
            if not tid or not gid:
                continue
            if feature == 'transcript':
                ref = is_reference(tid, gid, source, attrs)
                if (not ref) and passing and tid not in passing:
                    continue
                transcript_meta[tid] = (chrom, strand, gid, source, ref)
            elif feature == 'exon':
                exon_map[tid].append(Exon(int(start), int(end)))

    transcripts: List[Transcript] = []
    for tid, meta in transcript_meta.items():
        if tid not in exon_map:
            continue
        chrom, strand, gid, source, ref = meta
        transcripts.append(Transcript(chrom, strand, normalize_exons(exon_map[tid]), gid, tid, source, sample_id, ref))
    return transcripts


def junctions_close(a: List[Tuple[int, int]], b: List[Tuple[int, int]], tol: int) -> bool:
    if len(a) != len(b):
        return False
    return all(abs(x1 - x2) <= tol and abs(y1 - y2) <= tol for (x1, y1), (x2, y2) in zip(a, b))


def transcripts_match(a: Transcript, b: Transcript, end_tol: int, junction_tol: int) -> bool:
    if a.chrom != b.chrom or a.strand != b.strand or a.exon_count != b.exon_count:
        return False
    if not junctions_close(a.junctions(), b.junctions(), junction_tol):
        return False
    a5, a3 = a.ends()
    b5, b3 = b.ends()
    return abs(a5 - b5) <= end_tol and abs(a3 - b3) <= end_tol


def collapse_novel(novel_txs: List[Transcript], end_tol: int, junction_tol: int, min_samples: int) -> List[Transcript]:
    buckets: Dict[Tuple[str, str, int], List[Transcript]] = defaultdict(list)
    for tx in novel_txs:
        buckets[(tx.chrom, tx.strand, tx.exon_count)].append(tx)

    collapsed: List[Transcript] = []
    tx_counter = 1
    gene_counter = 1
    for _, bucket in buckets.items():
        groups: List[List[Transcript]] = []
        for tx in bucket:
            placed = False
            for group in groups:
                if transcripts_match(tx, group[0], end_tol, junction_tol):
                    group.append(tx)
                    placed = True
                    break
            if not placed:
                groups.append([tx])
        for group in groups:
            sample_count = len({t.sample_id for t in group})
            if sample_count < min_samples:
                continue
            rep = group[0]
            gene_id = f'gCOLL_{gene_counter:06d}'
            transcript_id = f'tCOLL_{tx_counter:06d}'
            tx_counter += 1
            gene_counter += 1
            collapsed.append(Transcript(rep.chrom, rep.strand, rep.exons, gene_id, transcript_id, 'IsoQuant', 'consensus', False))
    return collapsed


def write_gtf(transcripts: List[Transcript], out_gtf: str) -> None:
    with open(out_gtf, 'w') as handle:
        handle.write('# Consensus transcript annotation\n')
        gene_written: set[Tuple[str, str]] = set()
        for tx in sorted(transcripts, key=lambda t: (t.chrom, t.strand, t.start, t.end, t.transcript_id)):
            gene_key = (tx.gene_id, tx.strand)
            if gene_key not in gene_written:
                gene_written.add(gene_key)
                handle.write('\t'.join([tx.chrom, tx.source, 'gene', str(tx.start), str(tx.end), '.', tx.strand, '.', f'gene_id "{tx.gene_id}";']) + '\n')
            handle.write('\t'.join([tx.chrom, tx.source, 'transcript', str(tx.start), str(tx.end), '.', tx.strand, '.', f'gene_id "{tx.gene_id}"; transcript_id "{tx.transcript_id}";']) + '\n')
            for i, exon in enumerate(tx.exons, start=1):
                handle.write('\t'.join([tx.chrom, tx.source, 'exon', str(exon.start), str(exon.end), '.', tx.strand, '.', f'gene_id "{tx.gene_id}"; transcript_id "{tx.transcript_id}"; exon_number "{i}";']) + '\n')


def main() -> None:
    parser = argparse.ArgumentParser(description='Collapse per-sample transcript GTFs into a consensus GTF.')
    parser.add_argument('--samples', required=True)
    parser.add_argument('--project-dir', required=True)
    parser.add_argument('--output-gtf', required=True)
    parser.add_argument('--stats-output', required=True)
    parser.add_argument('--gtf-name', default='OUT.transcript_models.gtf')
    parser.add_argument('--matrix-name', default='transcript_matrix.csv')
    parser.add_argument('--min-reads', type=int, default=MIN_READS_THRESHOLD)
    parser.add_argument('--min-samples', type=int, default=MIN_SAMPLES_FOR_NOVEL)
    parser.add_argument('--end-tol', type=int, default=100)
    parser.add_argument('--junction-tol', type=int, default=10)
    args = parser.parse_args()

    samples = resolve_sample_dirs(args.project_dir, load_samples_tsv(args.samples))
    all_txs: List[Transcript] = []
    sample_stats: List[dict] = []
    for sample in samples:
        txs = read_gtf(
            str(Path(sample['out_dir']) / args.gtf_name),
            sample['sample_id'],
            args.min_reads,
            str(Path(sample['out_dir']) / args.matrix_name),
        )
        all_txs.extend(txs)
        sample_stats.append({'sample_id': sample['sample_id'], 'transcripts_after_read_filter': len(txs)})

    ref_txs = [t for t in all_txs if t.is_reference]
    novel_txs = [t for t in all_txs if not t.is_reference]
    ref_by_tid = {t.transcript_id: t for t in ref_txs}
    collapsed_novel = collapse_novel(novel_txs, args.end_tol, args.junction_tol, args.min_samples)
    final_txs = list(ref_by_tid.values()) + collapsed_novel

    ensure_dir(Path(args.output_gtf).parent)
    write_gtf(final_txs, args.output_gtf)

    with open(args.stats_output, 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['metric', 'value'])
        writer.writerow(['reference_transcripts', len(ref_by_tid)])
        writer.writerow(['novel_transcripts_before_collapse', len(novel_txs)])
        writer.writerow(['novel_transcripts_after_collapse', len(collapsed_novel)])
        writer.writerow(['total_transcripts_final', len(final_txs)])
        writer.writerow([])
        writer.writerow(['sample_id', 'transcripts_after_read_filter'])
        for row in sample_stats:
            writer.writerow([row['sample_id'], row['transcripts_after_read_filter']])


if __name__ == '__main__':
    main()
