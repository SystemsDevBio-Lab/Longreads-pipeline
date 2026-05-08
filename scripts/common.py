#!/usr/bin/env python3
"""Common utilities for the long-read scRNA-seq multi-sample pipeline."""

from __future__ import annotations

import csv
import gzip
import os
import re
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

from scipy import sparse
from scipy.io import mmwrite

BARCODE_UMI_RE = re.compile(r"^(?P<barcode>[A-Za-z0-9]{16}).*?_(?P<umi>[A-Za-z0-9]{12})#")


def ensure_dir(path: str | os.PathLike) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def open_maybe_gzip(path: str | os.PathLike, mode: str = "rt"):
    path = str(path)
    if path.endswith('.gz'):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_read_id(read_id: str) -> Tuple[str, str, str]:
    m = BARCODE_UMI_RE.search(read_id)
    if not m:
        raise ValueError(f"Could not parse barcode/UMI from read_id: {read_id}")
    barcode = m.group('barcode')
    umi = m.group('umi')
    return barcode, umi, barcode + umi


def load_samples_tsv(samples_tsv: str | os.PathLike) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    with open(samples_tsv, 'r', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'sample_id', 'fastq_path'}
        if not reader.fieldnames or not required.issubset(set(reader.fieldnames)):
            raise ValueError(
                f"samples.tsv must contain at least these columns: {sorted(required)}"
            )
        for row in reader:
            if not row.get('sample_id'):
                continue
            row = {k: (v.strip() if isinstance(v, str) else v) for k, v in row.items()}
            if not row.get('workdir'):
                row['workdir'] = ''
            rows.append(row)
    return rows


def resolve_sample_dirs(project_dir: str | os.PathLike, samples: List[Dict[str, str]]) -> List[Dict[str, str]]:
    project_dir = Path(project_dir)
    resolved = []
    for row in samples:
        sample_id = row['sample_id']
        workdir = (row.get('workdir') or '').strip()
        if not workdir:
            sample_dir = project_dir / 'per_sample' / sample_id
        else:
            workdir_path = Path(workdir)
            sample_dir = workdir_path if workdir_path.is_absolute() else project_dir / workdir_path
        row = dict(row)
        row['sample_dir'] = str(sample_dir)
        row['out_dir'] = str(sample_dir / 'OUT')
        resolved.append(row)
    return resolved


def write_text(path: str | os.PathLike, text: str) -> None:
    with open(path, 'w') as handle:
        handle.write(text)


def sorted_union(items: Iterable[str]) -> List[str]:
    return sorted({str(x) for x in items if str(x) != ''})


def write_lines(path: str | os.PathLike, values: Sequence[str]) -> None:
    with open(path, 'w') as handle:
        if values:
            handle.write("\n".join(str(v) for v in values) + "\n")
        else:
            handle.write("")


def reindex_matrix(matrix_df, row_names: Sequence[str] | None = None, col_names: Sequence[str] | None = None):
    result = matrix_df
    if row_names is not None:
        result = result.reindex(list(row_names), fill_value=0)
    if col_names is not None:
        result = result.reindex(columns=list(col_names), fill_value=0)
    if getattr(result, 'empty', False):
        return result
    return result.astype('int32')


def write_matrix_market_bundle(matrix, row_names: Sequence[str], col_names: Sequence[str], out_prefix: str | os.PathLike) -> None:
    out_prefix = Path(out_prefix)
    ensure_dir(out_prefix.parent)
    mmwrite(str(out_prefix.with_suffix('.mtx')), sparse.coo_matrix(matrix))
    write_lines(out_prefix.parent / f'{out_prefix.name}_rows.txt', list(row_names))
    write_lines(out_prefix.parent / f'{out_prefix.name}_cols.txt', list(col_names))
