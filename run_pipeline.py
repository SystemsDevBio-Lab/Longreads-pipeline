#!/usr/bin/env python3
"""One-command orchestrator for the long-read scRNA-seq multi-sample pipeline.

Changes in this version:
1. Sample-level parallel execution via --jobs.
2. Clean rerun of BLAZE to avoid stale intermediate files causing StopIteration.
3. Better logging with cwd/command written into each log file.
"""

from __future__ import annotations

import argparse
import csv
import shlex
import subprocess
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path

# Assumes scripts.common is available in the same project.
from scripts.common import ensure_dir, load_samples_tsv, resolve_sample_dirs


def run_command(command: str, log_path: Path, cwd: Path | None = None) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w") as log_handle:
        if cwd is not None:
            log_handle.write(f"[INFO] cwd: {cwd}\n")
        log_handle.write(f"[INFO] command: {command}\n")
        log_handle.flush()
        process = subprocess.run(
            ["bash", "-lc", command],
            cwd=str(cwd) if cwd is not None else None,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
        )
    if process.returncode != 0:
        raise RuntimeError(f"Command failed: {command}\nSee log: {log_path}")


def wrap_with_conda_env(command: str, env_name: str) -> str:
    escaped_command = command.replace("'", "'\\''")
    escaped_env = env_name.replace("'", "'\\''")
    return (
        "set -euo pipefail; "
        "if command -v conda >/dev/null 2>&1; then "
        "eval \"$(conda shell.bash hook)\"; "
        "elif [ -f \"$HOME/miniconda3/etc/profile.d/conda.sh\" ]; then "
        ". \"$HOME/miniconda3/etc/profile.d/conda.sh\"; "
        "elif [ -f \"$HOME/anaconda3/etc/profile.d/conda.sh\" ]; then "
        ". \"$HOME/anaconda3/etc/profile.d/conda.sh\"; "
        "else echo \"[ERR] conda not found\" >&2; exit 127; fi; "
        f"conda activate '{escaped_env}'; "
        f"{escaped_command}; "
        "conda deactivate"
    )


def run_parallel_samples(samples: list[dict[str, str]], jobs: int, worker_fn) -> None:
    if jobs <= 1:
        for sample in samples:
            worker_fn(sample)
        return

    with ThreadPoolExecutor(max_workers=jobs) as executor:
        future_to_sample = {
            executor.submit(worker_fn, sample): sample["sample_id"] for sample in samples
        }
        pending = set(future_to_sample)

        while pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                sample_id = future_to_sample[future]
                try:
                    future.result()
                    print(f"[OK] {sample_id} finished")
                except Exception:
                    for other in pending:
                        other.cancel()
                    raise RuntimeError(f"[FAILED] sample {sample_id}") from future.exception()


def remove_if_exists(path: Path) -> None:
    if path.exists():
        if path.is_file() or path.is_symlink():
            path.unlink()
        else:
            raise IsADirectoryError(f"Expected file but got directory: {path}")


def cleanup_blaze_outputs(sample_dir: Path, sample_id: str) -> None:
    """Remove stale BLAZE outputs that can cause mismatched intermediate state.

    The observed StopIteration in BLAZE happens after:
    - 'Search barcode in reads' is skipped because prior output exists
    - later step reads fewer barcode records than expected for current read batches

    This usually means old/incomplete intermediate files are being reused.
    """
    candidates = [
        sample_dir / f"{sample_id}_matched_reads.fastq.gz",
        sample_dir / f"{sample_id}putative_bc.csv",
        sample_dir / f"{sample_id}whitelist.csv",
        sample_dir / f"{sample_id}emtpy_bc_list.csv",  # BLAZE writes this misspelled name
        sample_dir / f"{sample_id}knee_plot.png",
        sample_dir / f"{sample_id}.barcode_count.tsv",
        sample_dir / f"{sample_id}.barcode_rank.tsv",
    ]
    for path in candidates:
        remove_if_exists(path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the full multi-sample long-read scRNA-seq pipeline."
    )
    parser.add_argument("--project-dir", required=True)
    parser.add_argument("--samples", required=True, help="samples.tsv")
    parser.add_argument("--reference", required=True, help="Reference genome FASTA")
    parser.add_argument(
        "--annotation", required=True, help="Reference GTF or IsoQuant gene DB path"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=16,
        help="Total CPU threads budget for the whole pipeline",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of samples to process in parallel",
    )
    parser.add_argument("--expect-cells", type=int, default=30000)
    parser.add_argument(
        "--blaze-conda-env", default="blaze_env", help="Conda environment name for blaze"
    )
    parser.add_argument(
        "--isoquant-conda-env",
        default="isoquant",
        help="Conda environment name for isoquant",
    )
    parser.add_argument(
        "--force-clean-blaze",
        action="store_true",
        default=True,
        help="Clean stale BLAZE intermediate outputs before rerunning BLAZE",
    )
    parser.add_argument("--skip-blaze", action="store_true")
    parser.add_argument("--skip-alignment", action="store_true")
    parser.add_argument("--skip-first-pass-isoquant", action="store_true")
    parser.add_argument("--skip-gene-matrix", action="store_true")
    parser.add_argument("--skip-bu", action="store_true")
    parser.add_argument("--skip-filter-merge", action="store_true")
    parser.add_argument("--skip-collapse", action="store_true")
    parser.add_argument("--skip-create-consensus-db", action="store_true")
    parser.add_argument("--skip-second-pass-isoquant", action="store_true")
    parser.add_argument("--skip-transcript-matrix", action="store_true")
    args = parser.parse_args()

    if args.threads < 1:
        raise ValueError("--threads must be >= 1")
    if args.jobs < 1:
        raise ValueError("--jobs must be >= 1")

    project_dir = Path(args.project_dir)
    samples = resolve_sample_dirs(project_dir, load_samples_tsv(args.samples))
    jobs = min(args.jobs, len(samples)) if samples else 1
    per_sample_threads = max(1, args.threads // jobs)

    logs_dir = ensure_dir(project_dir / "multisample" / "logs")
    matrices_dir = ensure_dir(project_dir / "multisample" / "matrices")
    merged_dir = ensure_dir(project_dir / "multisample" / "merged_bam")
    consensus_dir = ensure_dir(project_dir / "multisample" / "consensus_gtf")
    reannotation_dir = ensure_dir(project_dir / "multisample" / "reannotation")
    consensus_db_dir = ensure_dir(project_dir / "multisample" / "consensus_db")

    print(f"[INFO] total threads budget: {args.threads}")
    print(f"[INFO] parallel sample jobs: {jobs}")
    print(f"[INFO] threads per sample job: {per_sample_threads}")

    manifest_out = project_dir / "metadata" / "resolved_samples.tsv"
    ensure_dir(manifest_out.parent)
    with open(manifest_out, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample_id", "fastq_path", "workdir", "sample_dir", "out_dir"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(samples)

    def process_one_sample(sample: dict[str, str]) -> None:
        sample_id = sample["sample_id"]
        sample_dir = Path(sample["sample_dir"])
        out_dir = Path(sample["out_dir"])
        ensure_dir(sample_dir)
        ensure_dir(out_dir)

        matched_fastq = sample_dir / f"{sample_id}matched_reads.fastq.gz"
        sam_path = sample_dir / f"{sample_id}.sam"
        bam_path = sample_dir / f"{sample_id}.bam"
        sorted_bam = sample_dir / f"{sample_id}_sorted.bam"

        if not args.skip_blaze:
            if args.force_clean_blaze:
                cleanup_blaze_outputs(sample_dir, sample_id)

            blaze_cmd = (
                f"blaze --threads {per_sample_threads} --expect-cell {args.expect_cells} "
                f"--overwrite "
                f"--output-prefix {shlex.quote(sample_id)} "
                f"--output-fastq matched_reads.fastq.gz "
                f"{shlex.quote(sample['fastq_path'])}"
            )
            run_command(
                wrap_with_conda_env(blaze_cmd, args.blaze_conda_env),
                logs_dir / f"{sample_id}.blaze.log",
                cwd=sample_dir,
            )

        if not args.skip_alignment:
            cmd = (
                f"minimap2 -ax splice -t {per_sample_threads} -k14 --secondary=no "
                f"{shlex.quote(args.reference)} {shlex.quote(str(matched_fastq))} "
                f"> {shlex.quote(str(sam_path))}"
            )
            run_command(cmd, logs_dir / f"{sample_id}.minimap2.log", cwd=sample_dir)
            run_command(
                f"samtools view -@ {per_sample_threads} -bS {shlex.quote(str(sam_path))} > {shlex.quote(str(bam_path))}",
                logs_dir / f"{sample_id}.samtools_view.log",
                cwd=sample_dir,
            )
            run_command(
                f"samtools sort -@ {per_sample_threads} {shlex.quote(str(bam_path))} -o {shlex.quote(str(sorted_bam))}",
                logs_dir / f"{sample_id}.samtools_sort.log",
                cwd=sample_dir,
            )
            run_command(
                f"samtools index -@ {per_sample_threads} {shlex.quote(str(sorted_bam))}",
                logs_dir / f"{sample_id}.samtools_index.log",
                cwd=sample_dir,
            )
            remove_if_exists(sam_path)
            remove_if_exists(bam_path)

        if not args.skip_first_pass_isoquant:
            isoquant_cmd = (
                f"isoquant.py --bam {shlex.quote(str(sorted_bam))} "
                f"--reference {shlex.quote(args.reference)} "
                f"--genedb {shlex.quote(args.annotation)} "
                f"-t {per_sample_threads} --data_type nanopore "
                f"--model_construction_strategy default_ont "
                f"-o {shlex.quote(str(out_dir))}"
            )
            run_command(
                wrap_with_conda_env(isoquant_cmd, args.isoquant_conda_env),
                logs_dir / f"{sample_id}.isoquant_firstpass.log",
                cwd=sample_dir,
            )

        if not args.skip_bu:
            cmd = (
                f"python {shlex.quote(str(project_dir / 'scripts' / 'transcripts_bu.py'))} "
                f"--transcript-reads {shlex.quote(str(out_dir / 'OUT.transcript_model_reads.tsv.gz'))} "
                f"--bam {shlex.quote(str(sorted_bam))} "
                f"--output-dir {shlex.quote(str(out_dir))} "
                f"--sample-id {shlex.quote(sample_id)}"
            )
            run_command(cmd, logs_dir / f"{sample_id}.bu.log", cwd=sample_dir)

    run_parallel_samples(samples, jobs, process_one_sample)

    if not args.skip_gene_matrix:
        cmd = (
            f"python {shlex.quote(str(project_dir / 'scripts' / 'gene_quantification.py'))} "
            f"--samples {shlex.quote(str(manifest_out))} --project-dir {shlex.quote(str(project_dir))}"
        )
        run_command(cmd, logs_dir / "gene_quantification.log", cwd=project_dir)

    if not args.skip_filter_merge:
        cmd = (
            f"python {shlex.quote(str(project_dir / 'scripts' / 'filter_and_merge_bam.py'))} "
            f"--samples {shlex.quote(str(manifest_out))} "
            f"--project-dir {shlex.quote(str(project_dir))} "
            f"--output-dir {shlex.quote(str(merged_dir))} "
            f"--threads {args.threads}"
        )
        run_command(cmd, logs_dir / "filter_and_merge_bam.log", cwd=project_dir)

    consensus_gtf = consensus_dir / "integrated_output_filtered.gtf"
    consensus_stats = consensus_dir / "integrated_statistics.csv"
    if not args.skip_collapse:
        cmd = (
            f"python {shlex.quote(str(project_dir / 'scripts' / 'collapse_gtf.py'))} "
            f"--samples {shlex.quote(str(manifest_out))} "
            f"--project-dir {shlex.quote(str(project_dir))} "
            f"--output-gtf {shlex.quote(str(consensus_gtf))} "
            f"--stats-output {shlex.quote(str(consensus_stats))}"
        )
        run_command(cmd, logs_dir / "collapse_gtf.log", cwd=project_dir)

    consensus_db = consensus_db_dir / "integrated_output_filtered.db"
    if not args.skip_create_consensus_db:
        cmd = (
            f"python {shlex.quote(str(project_dir / 'scripts' / 'create_db.py'))} "
            f"{shlex.quote(str(consensus_gtf))} "
            f"--db {shlex.quote(str(consensus_db))} --force"
        )
        run_command(cmd, logs_dir / "create_consensus_db.log", cwd=project_dir)

    merged_bam = merged_dir / "all_samples_merged.bam"
    second_pass_out = reannotation_dir / "OUT"
    if not args.skip_second_pass_isoquant:
        isoquant_cmd = (
            f"isoquant.py --bam {shlex.quote(str(merged_bam))} "
            f"--reference {shlex.quote(args.reference)} "
            f"--genedb {shlex.quote(str(consensus_db))} "
            f"-t {args.threads} --data_type nanopore "
            f"--model_construction_strategy default_ont "
            f"-o {shlex.quote(str(second_pass_out))}"
        )
        run_command(
            wrap_with_conda_env(isoquant_cmd, args.isoquant_conda_env),
            logs_dir / "isoquant_secondpass.log",
            cwd=project_dir,
        )

    if not args.skip_transcript_matrix:
        cmd = (
            f"python {shlex.quote(str(project_dir / 'scripts' / 'transcript_quantification.py'))} "
            f"--map-tsv {shlex.quote(str(second_pass_out / 'OUT.transcript_model_reads.tsv.gz'))} "
            f"--samples {shlex.quote(str(manifest_out))} "
            f"--project-dir {shlex.quote(str(project_dir))} "
            f"--outdir {shlex.quote(str(matrices_dir))} "
            f"--shared-barcodes {shlex.quote(str(matrices_dir / 'shared_barcodes.txt'))}"
        )
        run_command(cmd, logs_dir / "transcript_quantification.log", cwd=project_dir)

    print("Pipeline finished successfully.")
    print(f"Gene matrix (csv.gz): {matrices_dir / 'gene_matrix.csv.gz'}")
    print(f"Gene matrix (csv): {matrices_dir / 'gene_matrix.csv'}")
    print(f"Gene matrix (mtx): {matrices_dir / 'gene_matrix.mtx'}")
    print(f"Transcript matrix (npz): {matrices_dir / 'matrix.npz'}")
    print(f"Transcript matrix (csv.gz): {matrices_dir / 'transcript_matrix.csv.gz'}")
    print(f"Transcript matrix (csv): {matrices_dir / 'transcript_matrix.csv'}")
    print(f"Transcript matrix (mtx): {matrices_dir / 'transcript_matrix.mtx'}")
    print(f"Shared barcodes: {matrices_dir / 'shared_barcodes.txt'}")
    print(f"Consensus GTF: {consensus_gtf}")
    print(f"Consensus DB: {consensus_db}")


if __name__ == "__main__":
    main()
