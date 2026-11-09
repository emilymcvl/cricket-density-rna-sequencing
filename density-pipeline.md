# RNA-Seq Analysis Pipeline Log

**Author:** Emily McColville  
**Date:** May 26, 2026  
**Project:** Effects of Crowding Density on Gene Expression in *Acheta domesticus*  
**Description:** Differential gene expression analysis of Illumina RNA-seq data from low- and high-density reared crickets.

This log documents all terminal commands and R scripts used for data processing, alignment, differential expression analysis, and functional annotation. All commands are provided for reproducibility.


---

## Table of Contents

- [RNA-Seq Analysis Pipeline Log](#rna-seq-analysis-pipeline-log)
  - [Table of Contents](#table-of-contents)
  - [1. Environment Setup](#1-environment-setup)
  - [2. Project Setup \& Directory Structure](#2-project-setup--directory-structure)
  - [3. Organizing Data](#3-organizing-data)
  - [4. Quality Control and Trimming](#4-quality-control-and-trimming)
    - [Pre-trimming QC](#pre-trimming-qc)
    - [Trimming Adapters](#trimming-adapters)
    - [Post-trimming QC](#post-trimming-qc)
    - [Fix Sample Naming](#fix-sample-naming)
  - [5. Sample Design File](#5-sample-design-file)
  - [6. Aligning to Transcriptome](#6-aligning-to-transcriptome)
    - [Download and Prepare Reference](#download-and-prepare-reference)
    - [Build Salmon Index](#build-salmon-index)
  - [7. Quantification with Salmon](#7-quantification-with-salmon)
    - [Check Mapping Rates](#check-mapping-rates)
  - [8. Combining Counts](#8-combining-counts)
  - [9. Differential Expression Analysis](#9-differential-expression-analysis)
    - [Design and model](#design-and-model)
    - [What the script does](#what-the-script-does)
    - [Output](#output)
    - [Output Columns](#output-columns)
  - [10. Functional Annotation with eggNOG-mapper](#10-functional-annotation-with-eggnog-mapper)
    - [Run eggNOG-mapper](#run-eggnog-mapper)
  - [11. Joining Annotations to DE Results](#11-joining-annotations-to-de-results)
    - [Output Files](#output-files)
  - [12. Directional GO \& KEGG Enrichment](#12-directional-go--kegg-enrichment)
    - [What the script does](#what-the-script-does-1)
    - [Output Files](#output-files-1)
  - [13. Annotating GO Terms with Descriptions](#13-annotating-go-terms-with-descriptions)
    - [Output Files](#output-files-2)
  - [14. Reducing GO Redundancy with rrvgo](#14-reducing-go-redundancy-with-rrvgo)
    - [What the script does](#what-the-script-does-2)
    - [Output Files](#output-files-3)
  - [15. GO Term Synthesis Summary](#15-go-term-synthesis-summary)
    - [Scope note](#scope-note)
    - [Output File](#output-file)
  - [16. Running the Full Downstream Pipeline](#16-running-the-full-downstream-pipeline)
  - [Appendix](#appendix)
    - [A1. Column Meanings (Full Reference)](#a1-column-meanings-full-reference)
    - [A2. eggNOG-mapper vs. BLAST](#a2-eggnog-mapper-vs-blast)
    - [A3. Manual BLAST (If Needed)](#a3-manual-blast-if-needed)
    - [A4. Notes on Annotation Strategy](#a4-notes-on-annotation-strategy)
    - [A5. Downstream Output Directory Tree](#a5-downstream-output-directory-tree)

---

## 1. Environment Setup

We use the following conda environments for this analysis:

| Environment | Purpose |
|-------------|---------|
| `bioinfo` | General bioinformatics tools |
| `fastqc` | FastQC and MultiQC quality control |
| `trim` | Trim Galore adapter trimming |
| `salmon` | Salmon quantification |
| `stats` | R scripts for DE analysis and visualization |
| `eggnog` | eggNOG-mapper functional annotation |
| `trinotate` | Transcript to protein prediction (TransDecoder) |

---

## 2. Project Setup & Directory Structure

Create the main project directory and subdirectories:

```bash
mkdir -p ~/work/density/{compressed-data,raw-data,merged-data,trimmed-data,refs,code,fastqc,clean-fastqc,bam,salmon,idx,eggnog}
cd ~/work/density
```

**Directory structure:**

```
~/work/density/
├── compressed-data/     # Original .ora files
├── raw-data/            # Decompressed .fastq.gz (all original lanes)
├── merged-data/         # Merged multi-lane + single-lane .fastq.gz (before trimming)
├── trimmed-data/        # Cleaned data .fastq.gz (after trimming)
├── bam/                 # Sorted, indexed BAM files (alignment output)
├── refs/                # genome.fa + transcripts.fa + predicted proteins
├── salmon/              # Salmon quantification results (per-sample)
├── idx/                 # Salmon index
├── code/                # Scripts, design.csv, sample IDs
├── fastqc/              # FastQC reports (pre-trimming)
├── clean-fastqc/        # FastQC reports (post-trimming)
└── eggnog/              # Annotation and enrichment analysis
```

---

## 3. Organizing Data

Decompress all `.ora` files to `.fastq.gz`:

```bash
cd ~/work/density
orad --gz compressed-data/*.fastq.ora
mv compressed-data/*.fastq.gz raw-data/
```

Merge multi-lane samples (E2 and E4) and move to `merged-data/`. Copy single-lane samples (A4, A5, A8, A10, E9, E10) to `merged-data/`:

```bash
cd ~/work/density

# Merge E2 (example)
zcat raw-data/E2_*_R1*.fastq.gz | gzip > merged-data/E2_R1.fastq.gz
zcat raw-data/E2_*_R2*.fastq.gz | gzip > merged-data/E2_R2.fastq.gz

# Copy single-lane samples
cp raw-data/A*_R*.fastq.gz merged-data/
```

---

## 4. Quality Control and Trimming

### Pre-trimming QC

Run FastQC on merged data:

```bash
conda activate fastqc
cd ~/work/density

fastqc merged-data/*.fastq.gz -o fastqc/
multiqc fastqc/ -o fastqc/
```

### Trimming Adapters

Trim adapters using Trim Galore and move trimmed data to `clean-data/`:

```bash
conda activate trim
cd ~/work/density/clean-data

trim_galore \
  --paired \
  --illumina \
  --length 50 \
  --quality 20 \
  --gzip \
  ../merged-data/*.fastq.gz
```

### Post-trimming QC

Run FastQC again on cleaned data:

```bash
conda activate fastqc
cd ~/work/density

fastqc clean-data/*.fastq.gz -o clean-fastqc/
multiqc clean-fastqc/ -o clean-fastqc/
```

### Fix Sample Naming

Standardize sample names for downstream processing:

```bash
cd ~/work/density/clean-data

# Example: A4_L_DEN_S77_L007_R1_001.fastq.gz → A4_R1.fastq.gz
# (Run rename script or manually via cp)
```

---

## 5. Sample Design File

Create `design.csv` in the `/code` directory to track sample names and conditions:

```csv
sample,condition
A10,low
A5,low
A4,low
A8,low
E10,high
E9,high
E2,high
E4,high
```

Extract sample IDs into `ids.txt` (no header):

```
A10
A5
A4
A8
E10
E9
E2
E4
```

This file is used to automate downstream analyses.

---

## 6. Aligning to Transcriptome


### Download and Prepare Reference

Download the *Acheta domesticus* reference transcriptome from [NCBI (GHVX00000000)](https://www.ncbi.nlm.nih.gov/Traces/wgs/?display=download). Rename and place in `refs/`:

```bash
cd ~/work/density/refs
# Download transcripts.fa manually (NCBI link above)
```

### Build Salmon Index

```bash
REF=refs/transcripts.fa
IDX=idx/salmon.idx

salmon index -t ${REF} -i ${IDX}
```

---

## 7. Quantification with Salmon

Run Salmon classification-based quantification for all samples. This automates the process using sample IDs from `ids.txt`:

```bash
conda activate salmon
cd ~/work/density

mkdir -p salmon

cat code/ids.txt | parallel -j 4 \
  "salmon quant \
    -i ${IDX} \
    -l A \
    --validateMappings \
    -1 clean-data/{}_R1.fq.gz \
    -2 clean-data/{}_R2.fq.gz \
    -o salmon/{}"
```

### Check Mapping Rates

```bash
grep -H "Mapping rate" salmon/*/logs/salmon_quant.log | \
  sed 's|.*/salmon/\([^/]*\)/logs/salmon_quant\.log:.*Mapping rate = \([^%]*\)%.*|salmon/\1 Mapping rate = \2%|'
```

Expected output (all samples should be >90%):

```
salmon/A10 Mapping rate = 96.5%
salmon/A4  Mapping rate = 95.3%
salmon/A5  Mapping rate = 96.3%
salmon/A8  Mapping rate = 96.6%
salmon/E10 Mapping rate = 95.5%
salmon/E2  Mapping rate = 92.1%
salmon/E4  Mapping rate = 91.6%
salmon/E9  Mapping rate = 94.0%
```

---

## 8. Combining Counts

Combine all per-sample quantification files into a single count matrix:

```bash
conda activate stats
cd ~/work/density

Rscript code/combine_transcripts.r
```

This produces `counts.csv` — a transcript × sample count matrix with all samples.

---

## 9. Differential Expression Analysis

Run DESeq2 to identify differentially expressed genes between low- and high-density conditions:

```bash
conda activate stats
cd ~/work/density/code

Rscript deseq2_density.R
```

### Design and model

The density experiment has a single experimental factor, `condition`, with two levels. **`low` is set as the reference level**, so the model is:

```
~ condition
```

and the single extracted coefficient is `condition_high_vs_low`. A positive `log2FoldChange` therefore means a gene is **up in high density** relative to low density.

### What the script does

1. Reads `design.csv` and the count matrix `../counts.csv`, rounding counts to integers.
2. Applies a low-count pre-filter (`filter_low_counts`): a transcript is kept only if it has ≥ 10 counts in at least as many samples as the smallest group (here, 4). This removes noise-dominated transcripts before testing.
3. Fits the DESeq2 model and extracts the `condition_high_vs_low` coefficient.
4. Applies apeglm LFC shrinkage to the coefficient. Shrinkage improves the log2FC estimate while the hypothesis test and FDR come from the standard Wald test for the same coefficient.
5. Formats and writes the results.

### Output

This produces `deseq2_results_density/results_high_vs_low.csv`.

### Output Columns

| Column | Description | Notes |
|--------|-------------|-------|
| `name` | Transcript ID (e.g., GHUU01000002.1) | Row identifier |
| `baseMean` | Average normalized expression (all samples) | Overall expression level |
| `foldChange` | Linear fold change (2^log2FC); high/low | Easy interpretation |
| `log2FoldChange` | log₂(high/low); +ve = up in high, −ve = up in low | apeglm-shrunk effect size |
| `lfcSE` | Standard error of log2FoldChange | Estimate precision |
| `stat` | Wald test statistic | Test statistic |
| `PValue` | Raw (unadjusted) p-value | Not corrected for multiple tests |
| `PAdj` | p-value adjusted by Hochberg method | Alternative FDR estimate |
| `FDR` | False Discovery Rate (Benjamini–Hochberg) | Multiple-test correction |
| `falsePos` | Expected false positives up to this row (index × FDR) | FDR interpretation |
| `A10, A5, A4, A8` | Normalized counts per sample (low density) | Heatmap/expression data |
| `E10, E9, E2, E4` | Normalized counts per sample (high density) | Heatmap/expression data |

**Key thresholds:**

- **FDR < 0.05:** ~5% false discovery rate among called genes.
- **\|log2FC\| > 1:** At least 2-fold change (log₂(2) = 1; log₂(0.5) = −1).

---

## 10. Functional Annotation with eggNOG-mapper

eggNOG-mapper requires protein sequences. Convert transcripts to predicted proteins using TransDecoder:

```bash
cd ~/work/density/refs

TransDecoder.LongOrfs -t transcripts.fa
TransDecoder.Predict -t transcripts.fa
```

This produces `transcripts.fa.transdecoder.pep`.

### Run eggNOG-mapper

```bash
conda activate eggnog
cd ~/work/density/eggnog

export EGGNOG_DATA_DIR=$(pwd)
mkdir -p temp

emapper.py \
  -i ../refs/transcripts.fa.transdecoder.pep \
  -m diamond \
  --data_dir "${EGGNOG_DATA_DIR}" \
  --output acheta_transdecoder_eggnog \
  --output_dir ./ \
  --cpu 8 \
  --temp_dir ./temp
```

This produces `acheta_transdecoder_eggnog.emapper.annotations`.

---

## 11. Joining Annotations to DE Results

Join eggNOG annotations with DESeq2 results:

```bash
conda activate stats
cd ~/work/density/code

Rscript annotate_de_density.R
```

The script reads the eggNOG annotation table (`../eggnog/acheta_transdecoder_eggnog.emapper.annotations`), strips the TransDecoder `.p#` suffix from query IDs so they match transcript IDs, and left-joins annotations onto every DESeq2 result file it finds in `deseq2_results_density/` (auto-discovered by the `results_*.csv` pattern). Columns are reordered as DESeq2 stats → key annotation fields → everything else.

### Output Files

Written under `deseq2_results_density/annotated/`:

- **`full/annotated_full_high_vs_low.csv`:** All genes with DESeq2 stats + eggNOG annotations (left join, keeping all DE results).
- **`significant/annotated_sig_high_vs_low.csv`:** Significantly DE genes filtered by **FDR < 0.05** and **|log2FoldChange| > 1**.
- **`annotation_summary.csv`:** Per-contrast counts of total genes, annotated genes, and significant DE genes.

**Filtering rationale:**

- **FDR < 0.05** controls false discovery rate at 5%.
- **|log2FC| > 1** ensures at least 2-fold change in either direction.

---

## 12. Directional GO & KEGG Enrichment

Test whether biological functions (GO terms) or pathways (KEGG KOs) are over-represented among DE genes, **separately for up- and down-regulated genes**:

```bash
conda activate stats
cd ~/work/density/code

Rscript enrich_density.R
```

### What the script does

1. Loads DE genes from `annotated/significant/annotated_sig_high_vs_low.csv` (filtered: FDR < 0.05, |log2FC| > 1).
2. Builds the statistical background universe from all annotated genes in `annotated/full/` (any gene with a `Preferred_name`, `GOs`, or `KEGG_ko`).
3. Splits DE genes into UP (log2FC > 0; higher in high density) and DOWN (log2FC < 0; higher in low density).
4. Runs `clusterProfiler::enricher` for each direction:
   - GO enrichment runs only if ≥ 5 annotated genes are present in that direction (`minGSSize = 5`, `maxGSSize = 500`); a `Count >= 5` safety filter is applied to results.
   - KEGG KO enrichment runs if ≥ 3 annotated genes are present (`minGSSize = 3`, `maxGSSize = 300`).
   - Cutoffs: p < 0.05, q < 0.2.

### Output Files

Written under `deseq2_results_density/comprehensive_go_analysis_density/`:

- **`GO_UP/go_UP_high_vs_low.csv`** and **`GO_DOWN/go_DOWN_high_vs_low.csv`:** Enriched GO terms per direction.
- **`GO_UP_annotated/`** and **`GO_DOWN_annotated/`:** Same terms with the `GO_Name` column (consumed by the next step).
- **`KEGG_UP/`** and **`KEGG_DOWN/`:** Enriched KEGG KOs per direction.
- **`directional_enrichment_summary.csv`:** UP/DOWN gene counts and enriched-term counts.

**Note:** This performs enrichment of individual KEGG KOs, not full pathways (limited by annotation completeness for *A. domesticus*).

---

## 13. Annotating GO Terms with Descriptions


```bash
conda activate stats
cd ~/work/density/code

Rscript annotate_go_terms_density.R
```

For every GO ID in the `GO_UP/` and `GO_DOWN/` enrichment files, the script looks up:

- **`GO_Name`** : the human-readable term name (e.g., "immune response").
- **`GO_Category`** : Biological Process, Molecular Function, or Cellular Component.
- **`GO_Definition`**: the full definition of the term.

### Output Files

- **`GO_UP_annotated/go_annotated_high_vs_low.csv`** and **`GO_DOWN_annotated/go_annotated_high_vs_low.csv`:** Enrichment results with the three description columns prepended.
- **`go_annotation_summary.csv`:** Per-file term counts broken down by BP / MF / CC.

---

## 14. Reducing GO Redundancy with rrvgo

GO enrichment often returns many semantically overlapping terms. `rrvgo` clusters similar terms and keeps the most significant representative of each cluster, making the results far easier to read and plot:

```bash
conda activate stats
cd ~/work/density/code

Rscript rrvgo_reduce_density.R
```

### What the script does

- Processes the comprehensive **UP** and **DOWN** annotated GO directories.
- For each ontology (BP / MF / CC) separately, computes a semantic **similarity matrix** (`calculateSimMatrix`, `method = "Rel"`) and reduces it (`reduceSimMatrix`) at a similarity **threshold of 0.8**, scoring terms by `−log10(p.adjust)`.


### Output Files

- **`GO_UP_reduced/`** and **`GO_DOWN_reduced/`:** Non-redundant GO terms, each with a `parentTerm` and `cluster_size` (how many raw terms collapsed into that representative).
- **`rrvgo_reduction_summary.csv`:** Per-file input vs output term counts and percent reduction.

---



## 15. GO Term Synthesis Summary

Produce a single human-readable summary of the contrast that combines the raw and rrvgo-reduced GO terms in both directions:

```bash
conda activate stats
cd ~/work/density/code

Rscript summaries-for-your-convenience.R
```

### Scope note

The flight pipeline performed *cross-contrast* synthesis — a master GO catalog, per-term "biographies", gene trajectories, and a presence/absence heatmap — all of which exist to find patterns **across many contrasts of equal statistical power**. The density experiment has a **single contrast**, so those comparisons would trivially report "present in 1 of 1" and are omitted. What remains is the per-contrast summary writer, which is genuinely useful for one contrast.

### Output File

- **`summaries-for-your-convenience/01_individual_comparisons/high_vs_low_summary.txt`:** A formatted summary with the biological question, an interpretation guide, UP/DOWN gene counts, and both raw and rrvgo-reduced GO term blocks (the reduced block shows the `ClusterSz` for each representative).

---

## 16. Running the Full Downstream Pipeline

Once `counts.csv` and the eggNOG annotations exist, the entire downstream analysis runs in order:

```bash
conda activate stats
cd ~/work/density/code

Rscript 1_deseq2l_density.R        # 1. DE testing (high vs low)
Rscript 2_annotate_de_density.R   # 2. Join eggNOG annotations
Rscript 3_enrich_density.R        # 3. Directional GO & KEGG enrichment
Rscript 4_annotate_go_terms_density.R            # 4. Add GO names/definitions
Rscript 5_rrvgo_reduce_density.R                 # 5. Reduce GO redundancy
Rscript 6_summaries-for-your-convenience.R     # 6. Summary text file
```

Each step depends on the output of the previous one, so the order matters.



## Appendix

### A1. Column Meanings (Full Reference)

| Column | Origin | Description | Unit/Type | Use |
|--------|--------|-------------|-----------|-----|
| `name` | DESeq2 | Transcript ID | Ensembl-style | Row identifier |
| `baseMean` | DESeq2 | Avg normalized expression (all) | Norm. counts | Overall level |
| `log2FoldChange` | DESeq2 | log₂(high / low), apeglm-shrunk | log₂ scale | Effect size |
| `foldChange` | DESeq2 | 2^log2FC = high/low | Linear | Easy reading |
| `lfcSE` | DESeq2 | SE of log2FC | log₂ scale | Precision |
| `stat` | DESeq2 | Wald test stat | — | Test result |
| `PValue` | DESeq2 | Raw p-value | 0–1 | Uncorrected |
| `PAdj` | DESeq2 | Hochberg-corrected p | Scientific | Alt. correction |
| `FDR` | DESeq2 | False Discovery Rate | 0–1 | Multiple testing |
| `falsePos` | DESeq2 | Expected FP count | Integer | FDR interpretation |
| `Preferred_name` | eggNOG | Best-matching gene name | Text | Primary ID |
| `Description` | eggNOG | Functional description | Text | Biology |
| `GOs` | eggNOG | Gene Ontology terms | Semicolon-sep. | GO enrichment |
| `KEGG_ko` | eggNOG | KEGG Orthology ID | e.g., K00407 | Pathway mapping |
| `KEGG_Pathway` | eggNOG | Pathway name(s) | Text | Pathway enrichment |
| `COG_category` | eggNOG | Broad functional class | Letter code | Function group |
| `A10, A5, A4, A8` | DESeq2 | Normalized counts | Decimal | Expression matrix |
| `E10, E9, E2, E4` | DESeq2 | Normalized counts | Decimal | Expression matrix |

### A2. eggNOG-mapper vs. BLAST

Both search for homologous proteins, but with key differences:

| Feature | BLAST (blastx + manual GO lookup) | eggNOG-mapper (DIAMOND mode) |
|---------|-----------------------------------|------------------------------|
| **Core method** | Any homology match | Homology refined via orthology (OGs) |
| **Accuracy** | Lower (includes paralogous hits) | Higher (orthology-based filtering) |
| **Output** | Raw hits (E-value, %ID, coverage) | Curated: OGs, GO, KEGG, COGs, names |
| **Annotation transfer** | From any hit | From orthologs only |
| **Benchmarks** | — | 11–15% fewer false-positive GO terms vs. BLAST |

