# Transcriptomics raw data

Raw RNA-seq count matrices and sample metadata for the three studies analyzed in this paper.
ASOs are identified by the anonymized scheme used throughout the repository (see `analytic/treatments.tsv`);
proprietary numeric compound identifiers have been remapped out of the metadata.

```
transcriptomics/
├── 21-01ac/   21-01ac_counts.csv.gz      transcript-level counts (samples = columns)
│              21-01ac_metadata.csv        per-sample metadata (treatment = anonymized ASO)
├── 23-0afa/   23-0afa_counts.csv.gz       gene-level counts
│              23-0afa_metadata.csv        per-sample metadata
│              23-0afa_samples.tsv         animal → test article
└── 24-0432/   24-0432_counts.csv.gz       gene-level counts
               24-0432_metadata.csv        per-sample metadata
               24-0432_inoculum.tsv        sample → prion inoculum (RML/CBH)
               24-0432_variants.tsv        ASO variant/category metadata
               analysis.R                  provenance: DESeq2/GSEA pipeline (see header)
               output/                      created when analysis.R is run
```

Notes:
- Count-matrix column headers are sequencing sample IDs; gene/transcript identifiers are standard
  Ensembl/MGI symbols.
- The DESeq2 differential-expression results derived from these counts are in `analytic/`
  (`deseq2_contrasts.tsv.gz`, `all_contrasts_0432.tsv.gz`), which the figure script reads.
- Metadata are the original sequencing-core/LIMS exports, copied verbatim except that the
  treatment/test-article identifiers have been replaced with the anonymized ASO labels.
