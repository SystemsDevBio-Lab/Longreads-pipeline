# Long-read scRNA-seq multi-sample assembly and quantification pipeline

This repository contains the multi-sample assembly and quantification pipeline used for single-cell long-read RNA-seq data. The workflow starts from per-sample FASTQ files, builds a consensus transcript set across samples, and produces synchronized gene- and transcript-level matrices.

## Overview

```text
FASTQ
  -> BLAZE demultiplexing
  -> minimap2 alignment
  -> samtools BAM processing
  -> first-pass IsoQuant
     -> gene quantification
     -> transcript model collection
  -> consensus transcript collapse
  -> consensus database construction
  -> second-pass IsoQuant re-annotation
  -> final gene and transcript matrices
```

The workflow keeps per-sample processing separate until transcript models and filtered BAM files are merged for multi-sample analysis.

## Requirements

Recommended software:

```bash
conda create -n longreads_env python=3.10 -y
conda activate longreads_env
conda install -c bioconda minimap2 samtools -y
pip install pysam pandas numpy scipy

conda create -n blaze_env python=3.10 -y
conda activate blaze_env
pip install blaze
conda deactivate

conda create -n isoquant python=3.10 -y
conda activate isoquant
pip install isoquant
conda deactivate
```

Required inputs:

- reference genome FASTA
- reference annotation GTF or IsoQuant-compatible gene database
- one FASTQ file per sample
- sample sheet with `sample_id`, `fastq_path`, and optional `workdir`

## Project layout

```text
project/
├── metadata/
│   ├── samples.template.tsv
│   └── resolved_samples.template.tsv
├── reference/
│   ├── genome.fa
│   └── genes.gtf
├── scripts/
│   ├── common.py
│   ├── gene_quantification.py
│   ├── transcripts_bu.py
│   ├── filter_and_merge_bam.py
│   ├── collapse_gtf.py
│   ├── create_db.py
│   └── transcript_quantification.py
├── per_sample/
│   ├── sample1/
│   ├── sample2/
│   └── ...
├── multisample/
│   ├── merged_bam/
│   ├── consensus_gtf/
│   ├── consensus_db/
│   ├── reannotation/
│   ├── matrices/
│   └── logs/
└── run_pipeline.py
```

## Sample sheet

Use `metadata/samples.template.tsv` to define input files:

```tsv
sample_id	fastq_path	workdir
```

`sample_id` is used throughout the pipeline for file names and matrix columns. If `workdir` is empty, outputs are written to `per_sample/{sample_id}/`. Relative `workdir` values are resolved under `--project-dir`.

## Core outputs per sample

```text
per_sample/{sample_id}/
├── {sample_id}_matched_reads.fastq.gz
├── {sample_id}_sorted.bam
└── OUT/
    ├── OUT.read_assignments.tsv.gz
    ├── OUT.transcript_model_reads.tsv.gz
    ├── OUT.transcript_models.gtf
    ├── gene_matrix.csv.gz
    ├── OUT.transcript_model_reads_with_length.tsv.gz
    └── {sample_id}_bu.csv
```

## Main steps

### 1. Reads deconstruction

BLAZE identifies cell barcodes from long-read single-cell RNA-seq reads and writes deconstructed reads for each sample.

Output:

```text
{sample_id}_matched_reads.fastq.gz
```

### 2. Reads mapping

Reads are aligned to the reference genome with minimap2 and converted to sorted, indexed BAM files with samtools.

Outputs:

```text
{sample_id}_sorted.bam
{sample_id}_sorted.bam.bai
```

### 3. Transcripts annotation

Each sample is annotated independently with IsoQuant using the sorted BAM and reference annotation.

Key outputs:

```text
OUT.read_assignments.tsv.gz
OUT.transcript_model_reads.tsv.gz
OUT.transcript_models.gtf
```

### 4. Gene quantification

Gene counts are generated from `OUT.read_assignments.tsv.gz` after barcode-UMI deduplication. A barcode-UMI is retained only when all supporting reads map to one gene. Per-sample matrices are merged by the union of genes and concatenated cell columns.

Final output:

```text
multisample/matrices/gene_matrix.csv.gz
```

### 5. BAM filter and merge

Reads with transcript assignment not equal to `*` are retained from each BAM. Filtered BAMs are then merged across samples.

Outputs:

```text
multisample/merged_bam/
├── {sample_id}_filtered.bam
├── {sample_id}_filtered.bam.bai
├── filtering_statistics.tsv
├── all_samples_merged.bam
└── all_samples_merged.bam.bai
```

### 6. GTF filter and collapse

First-pass transcript GTF files are collapsed into a consensus transcript set. 

Outputs:

```text
multisample/consensus_gtf/
├── integrated_output_filtered.gtf
└── integrated_statistics.csv
```

### 7. Transcripts re-annotation

The merged BAM is re-annotated against the collapsed GTF.

Inputs:

```text
multisample/merged_bam/all_samples_merged.bam
multisample/consensus_db/integrated_output_filtered.db
```

Key output:

```text
multisample/reannotation/OUT/OUT.transcript_model_reads.tsv.gz
```

### 8. Transcripts quantification

The final transcript matrix is built from the transcripts re-annotation table and all per-sample barcode-UMI tables. 

Outputs:

```text
multisample/matrices/
├──transcript_matrix.csv.gz
├──transcript_matrix.csv
├──transcript_matrix.mtx
├──ranscript_matrix_rows.txt
├──transcript_matrix_cols.txt
├──matrix.npz
├──transcripts.txt
└──barcodes.txt
```

## Run the pipeline

```bash
python run_pipeline.py \
  --project-dir /path/to/project \
  --samples /path/to/project/metadata/samples.template.tsv \
  --reference /path/to/project/reference/genome.fa \
  --annotation /path/to/project/reference/genes.gtf \
  --threads 16
```

## Final deliverables

The main outputs are written under `multisample/matrices/` and `multisample/consensus_gtf/`:

```text
gene_matrix.csv.gz
transcript_matrix.csv.gz
transcript_matrix.mtx
transcript_matrix_rows.txt
transcript_matrix_cols.txt
shared_barcodes.txt
integrated_output_filtered.gtf
integrated_statistics.csv
```

Gene and transcript matrices use the same `sampleid_barcode` columns and the same column order. Missing entries are filled with zero.

