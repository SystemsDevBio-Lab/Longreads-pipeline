# Long-read scRNA-seq Multi-sample Pipeline

This repository contains a fully reorganized long-read single-cell RNA-seq workflow for Oxford Nanopore data. It keeps the original biological logic of your pipeline, but makes the implementation internally consistent, English-only, sample-sheet driven, and easier to run across all 25 samples.

## 1. Design goals

This revision enforces three principles across the entire workflow:

1. **All steps are sample-sheet driven.** No Python script contains hard-coded sample paths or a hard-coded 25-sample dictionary.
2. **All intermediate files follow a single directory convention.** Per-sample outputs are written under `per_sample/{sample_id}/`, while pooled results are written under `multisample/`.
3. **The gene branch and transcript branch are explicitly synchronized.** The same `sample_id`, `sampleid_barcode` naming rule, and standardized filenames are used throughout all downstream steps.

## 2. Workflow overview

The end-to-end workflow is:

```text
Raw FASTQ
  -> BLAZE
  -> minimap2
  -> samtools
  -> first-pass IsoQuant
     -> gene branch: per-sample gene matrices -> merged gene matrix
     -> transcript branch:
          per-sample BU tables
          -> BAM filtering
          -> pooled BAM merge
          -> consensus GTF collapse
          -> GTF to gffutils DB
          -> second-pass IsoQuant
          -> final transcript matrix
```

## 3. Software

Recommended environment:

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

## 4. Required inputs

You need the following files before running the full workflow:

- reference genome FASTA
- reference annotation GTF or IsoQuant-compatible gene database
- one FASTQ file per sample
- `metadata/samples.template.tsv`

## 5. Standard project layout

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
│   ├── 241218_01_HT/
│   ├── 241218_01_VT/
│   ├── ...
│   └── 250530_01_VL/
├── multisample/
│   ├── merged_bam/
│   ├── consensus_gtf/
│   ├── consensus_db/
│   ├── reannotation/
│   ├── matrices/
│   └── logs/
└── run_pipeline.py
```

## 6. The 25-sample manifest

Use the following 25 sample IDs consistently everywhere:

```text
241218_01_HT
241218_01_VT
250117_01_Head
250117_02_Trunk
250117_01_VL
250206_01_Head
250206_02_Trunk
250206_01_VL
250411_01_HT
250411_01_VT
250414_01_HT
250414_01_VT
nao2
neizang1
neizang2
qugan1
qugan2
sizhi1
sizhi2
241220_01_Head
241220_01_Trunk
241220_01_VL
250530_01_Head
250530_01_Trunk
250530_01_VL
```

A ready-to-edit template with your current FASTQ paths is provided in `metadata/samples.template.tsv`.

The required columns are:

```tsv
sample_id	fastq_path	workdir
```

Notes:

- `sample_id` must match the 25-sample list exactly.
- `fastq_path` is the absolute or relative path to the sample FASTQ.
- `workdir` may be left empty. If empty, the pipeline will automatically use `project/per_sample/{sample_id}`. Relative `workdir` values are resolved under `--project-dir`; absolute paths are kept as-is.

## 7. File naming conventions used by every script

For each sample, the pipeline expects or generates the following core files:

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

This naming scheme is shared across the documentation and all Python scripts.

## 8. Step-by-step logic

### 8.1 BLAZE

Each sample FASTQ is processed independently.

Output:

- `{sample_id}_matched_reads.fastq.gz`

### 8.2 minimap2 + samtools

Each matched FASTQ is aligned to the reference genome and converted to sorted/indexed BAM.

Output:

- `{sample_id}_sorted.bam`
- `{sample_id}_sorted.bam.bai`

### 8.3 First-pass IsoQuant

Each sample is processed independently using the sorted BAM and reference annotation.

Important first-pass outputs:

- `OUT.read_assignments.tsv.gz`
- `OUT.transcript_model_reads.tsv.gz`
- `OUT.transcript_models.gtf`

### 8.4 Gene quantification branch

The gene branch uses `OUT.read_assignments.tsv.gz` and applies barcode-UMI deduplication.

Rules:

1. Extract barcode and UMI from `read_id`.
2. Group reads by barcode-UMI.
3. Retain a barcode-UMI only if all supporting reads map to exactly one gene.
4. Count retained molecules into a per-sample gene-by-cell matrix.
5. Merge all per-sample matrices by taking the union of genes and concatenating columns.

Column naming:

```text
{sample_id}_{barcode}
```

Final pooled output:

```text
multisample/matrices/gene_matrix.csv.gz
```

### 8.5 Representative BU generation

For transcript-level analysis, each sample generates a representative BU table from:

- `OUT.transcript_model_reads.tsv.gz`
- `{sample_id}_sorted.bam`

Rules:

1. Keep reads with transcript assignment not equal to `*`.
2. Add `read_length` from BAM.
3. Group reads by barcode-UMI.
4. If one barcode-UMI maps to exactly one transcript, keep the longest read as the representative.
5. If one barcode-UMI maps to multiple transcripts, mark it as ambiguous.

Output:

```text
per_sample/{sample_id}/OUT/{sample_id}_bu.csv
```

Columns:

```text
bu,barcode,umi,read_id,transcript_id,read_length,status
```

### 8.6 BAM filtering and pooled merge

For each sample, the BAM is filtered using `OUT.transcript_model_reads.tsv.gz`.
Only reads with transcript assignment not equal to `*` are retained.

Outputs:

```text
multisample/merged_bam/
├── {sample_id}_filtered.bam
├── {sample_id}_filtered.bam.bai
├── filtering_statistics.tsv
├── all_samples_merged.bam
└── all_samples_merged.bam.bai
```

### 8.7 Consensus GTF construction

All first-pass `OUT.transcript_models.gtf` files are read and merged into a consensus transcript set.

Current rules in the generalized script:

- reference transcripts are preserved by transcript ID
- non-reference transcripts can be filtered by per-sample transcript matrix support
- novel transcripts are collapsed by exon count, chromosome, strand, transcript ends, and splice-junction tolerance
- a novel transcript must be observed in at least 2 samples to enter the consensus set

Default collapse tolerances:

- transcript end tolerance: 100 bp
- internal splice junction tolerance: 10 bp

Outputs:

```text
multisample/consensus_gtf/
├── integrated_output_filtered.gtf
└── integrated_statistics.csv
```

### 8.8 Consensus DB creation

The collapsed consensus GTF is converted into a gffutils SQLite database before reannotation.

Input:

- `multisample/consensus_gtf/integrated_output_filtered.gtf`

Output:

```text
multisample/consensus_db/
└── integrated_output_filtered.db
```

### 8.9 Second-pass IsoQuant

The pooled filtered BAM is re-annotated against the consensus database.

Input:

- `multisample/merged_bam/all_samples_merged.bam`
- `multisample/consensus_db/integrated_output_filtered.db`

Output directory:

```text
multisample/reannotation/OUT/
```

Key output used downstream:

```text
multisample/reannotation/OUT/OUT.transcript_model_reads.tsv.gz
```

### 8.10 Final transcript quantification

The final transcript matrix is generated from:

- second-pass `OUT.transcript_model_reads.tsv.gz`
- all per-sample `{sample_id}_bu.csv` tables

Rules:

1. Keep only read IDs that occur exactly once in the second-pass mapping table.
2. Discard all multi-mapped transcript assignments.
3. Map each retained read back to a `sampleid_barcode` column via BU tables.
4. Build a sparse transcript-by-cell matrix.

Outputs:

```text
multisample/matrices/
├── matrix.npz
├── transcripts.txt
└── barcodes.txt
```

## 9. One-command execution

After you place:

- FASTQ files
- reference genome
- annotation file
- scripts
- `samples.template.tsv`

into the expected structure, you can run the full pipeline with a single command:

```bash
python run_pipeline.py \
  --project-dir /path/to/project \
  --samples /path/to/project/metadata/samples.template.tsv \
  --reference /path/to/project/reference/genome.fa \
  --annotation /path/to/project/reference/genes.gtf \
  --threads 16
```

## 10. Recommended execution order inside `run_pipeline.py`

The integrated runner executes the following stages in order:

1. BLAZE
2. minimap2
3. samtools view/sort/index
4. first-pass IsoQuant
5. per-sample BU generation
6. gene matrix generation
7. BAM filtering and pooled merge
8. consensus GTF collapse
9. consensus DB creation
10. second-pass IsoQuant
11. final transcript quantification

## 11. Final expected deliverables

After successful completion, the most important files are:

```text
multisample/matrices/gene_matrix.csv.gz
multisample/matrices/matrix.npz
multisample/matrices/transcripts.txt
multisample/matrices/barcodes.txt
multisample/consensus_gtf/integrated_output_filtered.gtf
multisample/consensus_gtf/integrated_statistics.csv
```

## 12. Practical notes

- All scripts are now English-only.
- All scripts are CLI-based and reusable across projects.
- All scripts use `samples.template.tsv` or explicit arguments instead of internal hard-coded paths.
- The current generalized implementation assumes the BLAZE read header format still contains a 16 bp barcode and 12 bp UMI in the same pattern used by your original pipeline.
- The second-pass IsoQuant step now uses a gffutils SQLite DB created from the collapsed consensus GTF via `scripts/create_db.py`.



## Final merged matrix exports

This revision writes synchronized gene and transcript outputs under `multisample/matrices/`.

Shared column definition:

```text
shared_barcodes.txt
```

Gene matrix outputs:

```text
gene_matrix.csv.gz
gene_matrix.csv
gene_matrix.mtx
gene_matrix_rows.txt
gene_matrix_cols.txt
```

Transcript matrix outputs:

```text
transcript_matrix.csv.gz
transcript_matrix.csv
transcript_matrix.mtx
transcript_matrix_rows.txt
transcript_matrix_cols.txt
```

Both matrices use the same `sampleid_barcode` columns and identical column order. The shared barcode list is built from the union of gene-branch and transcript-BU barcodes, and missing entries are filled with 0.
