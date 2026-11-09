# cricket-density-rna-sequencing

# Hello

If you've found your way to this repo, hopefully you've enjoyed our paper! My name is Emily, and I was the lead author on this project.


![We hope you're excited!](images/meme.jpg)

The following repo includes all of the documents you need to follow along with the analysis for our study on how population density affects gene expression in crickets. The file `density-pipeline.md shows` the full pipeline we used. The original file paths are preserved to reflect how the analysis was run in practice, but to save space and avoid unnecessary clutter, I’ve only uploaded the key files and results here.

The code folder contains all of the scripts used for the analysis, in the order they were run, and the results folder contains the outputs of those scripts.

## What this project is

This project uses RNA-seq data to answer the question: how does rearing density affect gene expression in crickets? The analysis compares transcript abundance across the two density groups and then follows up by identifying enriched functional categories such as Gene Ontology (GO) terms and KEGG pathways.

This repo includes:

-   sample metadata and experimental design
-   count-matrix processing and normalization
-   differential expression analysis with DESeq2
-   functional annotation of genes via eggNOG
-   directional GO and KEGG enrichment analysis
-   summary reports for your convenience :) 



## Repository structure

``` text
.
├── README.md
├── density-pipeline-may-26.md
├── LICENSE
├── code/
│   ├── design.csv
│   ├── ids.txt
│   ├── 1_deseq2_multifactorial_density.R
│   ├── 2_annotate_de_multifactorial_density.R
│   ├── 3_enrich_multifactorial_density.R
│   ├── 4_annotate_go_terms_density.R
│   ├── 5_rrvgo_reduce_density.R
│   └── 6_summaries-for-your-convenience.R
├── results/
│   ├── annotated_sig_high_vs_low_ns.csv
│   ├── directional_enrichment_summary.csv
│   ├── GO_DOWN_annotated/
│   ├── GO_UP_annotated/
│   ├── KEGG_DOWN/
│   ├── KEGG_UP/
│   └── summaries-for-your-convenience/
└── ...
```



## Main workflow

The pipeline is organized as a standard RNA-seq analysis workflow:

1.  Load sample metadata and count data
2.  Build a DESeq2 model for the density contrast
3.  Filter low-count genes and normalize expression values
4.  Run differential expression analysis
5.  Annotate DE genes with functional information from eggNOG
6.  Test for directional enrichment of GO and KEGG terms
7.  Reduce redundant GO terms for cleaner biological interpretation
8.  Save summary files and report-ready outputs

## Key scripts

-   `1_deseq2_multifactorial_density.R`\
    Runs the DESeq2 differential expression analysis for the high-vs-low density comparison.

-   `2_annotate_de_multifactorial_density.R`\
    Joins gene annotations from eggNOG to the DE results so each gene has functional context.

-   `3_enrich_multifactorial_density.R`\
    Performs GO and KEGG enrichment for up-regulated and down-regulated genes.

-   `4_annotate_go_terms_density.R`\
    Adds descriptive names to GO IDs in enrichment results.

-   `5_rrvgo_reduce_density.R`\
    Reduces redundant GO terms using rrvgo to keep only the most informative clusters.

-   `6_summaries-for-your-convenience.R`\
    Produces human-readable summaries of the analysis outputs.

## Outputs

The `results/` directory contains downstream files such as:

-   differentially expressed genes
-   annotated DE results
-   GO and KEGG enrichment tables
-   reduced GO-term summaries
-   convenience summary reports for interpretation

These outputs are intended to support biological interpretation of which pathways and processes are altered by high-density rearing.

## Notes

This repo is built around a single factorial density contrast and is focused on transcriptome-level interpretation.. The scripts are set up for reproducible analysis of this cricket density experiment.

\