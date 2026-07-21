# analysis.R — ASO toxicity via DESeq2-based DGE scoring (study 24-0432)
#
# PROVENANCE / DOCUMENTATION ONLY — not part of the figure pipeline.
# Documents how the study-24 derived datasets used by the manuscript
# (analytic/umap_data.tsv, analytic/dge_scores_24.tsv, analytic/volcano_data.tsv,
# analytic/all_contrasts_0432.tsv.gz) were produced from the raw counts + metadata
# in transcriptomics/24-0432/. Requires extra packages beyond the figure script
# (DESeq2, uwot, fgsea, msigdbr, and the lab's `ericplot`) and writes its outputs to
# transcriptomics/24-0432/output/. ASOs use the anonymized labels in analytic/treatments.tsv.
#
# Computes G_tox continuous score (sum of signed L2FCs for G_tox genes).
# Runnable from any directory inside the repo (CWD is reset to repo root below).

library(tidyverse)
library(janitor)
library(DESeq2)
library(uwot)
library(fgsea)
library(msigdbr)

if (!requireNamespace("ericplot", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes", repos = "https://cloud.r-project.org")
  remotes::install_github("vallabhminikel/ericplot")
}

library(ericplot)

# Resolve repo root so all relative paths work regardless of CWD.
# Prefer the script's own location (Rscript invocation), fall back to walking
# up from CWD to the .git directory (interactive / sourced use).
.repo_root <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  if (length(file_arg) == 1) {
    rprojroot::find_root(rprojroot::is_git_root,
                         path = dirname(normalizePath(file_arg)))
  } else {
    rprojroot::find_root(rprojroot::is_git_root)
  }
})
setwd(.repo_root)
cat("Repo root:", .repo_root, "\n")

# ==============================================================================
# Step 1: Load raw counts + metadata for study 24
# ==============================================================================

cat("=== Loading data ===\n")

# ---- Study 24 counts (transcript → gene aggregation) ----
counts_24_raw <- read_csv("transcriptomics/24-0432/24-0432_counts.csv.gz",
                          show_col_types = FALSE)
counts_24 <- counts_24_raw %>%
  filter(!is.na(symbol) & symbol != "") %>%
  select(-gid, -tid) %>%
  group_by(symbol) %>%
  summarise(across(everything(), \(x) sum(x, na.rm = TRUE)), .groups = "drop")
count_mat_24 <- counts_24 %>% column_to_rownames("symbol") %>% as.matrix()
storage.mode(count_mat_24) <- "integer"

# ---- Study 24 metadata ----
meta_24_raw <- read_csv("transcriptomics/24-0432/24-0432_metadata.csv",
                        show_col_types = FALSE)
# Public data identifies ASOs by anonymized 'aso' column; restore the internal name.
meta_24_info <- read_tsv("transcriptomics/24-0432/24-0432_variants.tsv", show_col_types = FALSE) %>%
  rename(ionis_name = aso)

# ASO name + group lookup
aso_lookup_24 <- meta_24_info %>%
  filter(ionis_name != "saline") %>%
  transmute(
    ionis_id = as.character(ionis_name),
    aso_name = case_when(
      detail == "parent"    ~ "ASO 2",
      detail == "ASO 1"     ~ "ASO 1",
      detail %in% c("full 2'MOE", "2'OMe", "MsPA") ~ paste0("ASO 2 ", detail),
      detail %in% LETTERS   ~ paste0("ASO ", detail),
      TRUE ~ detail
    ),
    group = group
  )

# Condition categories (anonymized ASO identifiers)
TOXIC_ASOS     <- c("ASOC", "ASOD", "ASOF", "ASOI")
ACTIVE_ASO2    <- "ASO2"
MODIFIED_ASOS  <- c("ASO2-MOE", "ASO2-OMe2", "ASO2-MsPA1")
TOLERATED_ASOS <- meta_24_info %>%
  filter(group == "tolerated ASOs") %>% pull(ionis_name) %>% as.character()

meta_24 <- meta_24_raw %>%
  mutate(sample_id = as.character(ngsid),
         ionis_id  = as.character(treatment)) %>%
  left_join(aso_lookup_24, by = "ionis_id") %>%
  mutate(
    aso_name = ifelse(ionis_id == "PBS", "PBS", aso_name),
    aso_class = case_when(
      ionis_id == "PBS"              ~ "PBS",
      ionis_id %in% TOXIC_ASOS      ~ "toxic",
      ionis_id == ACTIVE_ASO2       ~ "toxic",
      ionis_id %in% MODIFIED_ASOS   ~ "modified",
      ionis_id %in% TOLERATED_ASOS  ~ "tolerated",
      TRUE                           ~ "other"
    ),
    study = "24-0432"
  ) %>%
  filter(sample_id %in% colnames(count_mat_24)) %>%
  select(sample_id, ionis_id, aso_name, aso_class, study)

meta_24 <- meta_24[match(colnames(count_mat_24), meta_24$sample_id), ]
rownames(meta_24) <- meta_24$sample_id
stopifnot(all(colnames(count_mat_24) == meta_24$sample_id))

cat(sprintf("Study 24: %d genes x %d samples\n", nrow(count_mat_24), ncol(count_mat_24)))

count_mat <- count_mat_24
meta <- meta_24
print(table(meta$aso_class))

# ==============================================================================
# Step 2: Inoculum assignment (ground truth from manifest)
# ==============================================================================

cat("\n=== Inoculum assignment ===\n")

inoculum_24 <- read_tsv("transcriptomics/24-0432/24-0432_inoculum.tsv", show_col_types = FALSE) %>%
  mutate(ngsid = as.character(ngsid))
meta_24 <- meta_24 %>%
  left_join(inoculum_24 %>% select(ngsid, inoculum), by = c("sample_id" = "ngsid"))
stopifnot(!any(is.na(meta_24$inoculum)))

cat(sprintf("Study 24 inoculum: %d RML, %d CBH\n",
            sum(meta_24$inoculum == "RML"), sum(meta_24$inoculum == "CBH")))

# ==============================================================================
# Step 3: DESeq2 — inoculum + toxic (study 24 only) → G_tox
# ==============================================================================

cat("\n=== DESeq2: ~ inoculum + pooled_condition (study 24) ===\n")

deseq_ids <- meta_24$sample_id[
  meta_24$ionis_id %in% c(TOXIC_ASOS, ACTIVE_ASO2, "PBS")
]
deseq_meta <- meta_24[match(deseq_ids, meta_24$sample_id), ] %>%
  mutate(
    pooled_condition = factor(
      ifelse(aso_class == "PBS", "PBS", "TOXIC"),
      levels = c("PBS", "TOXIC")
    ),
    inoculum = factor(inoculum, levels = c("CBH", "RML"))
  )
deseq_counts <- count_mat_24[, deseq_ids]

dds <- DESeqDataSetFromMatrix(
  countData = deseq_counts,
  colData   = deseq_meta,
  design    = ~ inoculum * pooled_condition
)
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
dds  <- DESeq(dds)

cat("DESeq2 coefficients:\n")
print(resultsNames(dds))

# Toxic effect in RML = main effect + interaction
res <- results(dds, contrast = list(
  c("pooled_condition_TOXIC_vs_PBS", "inoculumRML.pooled_conditionTOXIC")
))
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  filter(!is.na(padj))

# G_tox: genes DE in toxic vs PBS
g_tox <- res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 1, baseMean >= 10) %>%
  mutate(sign = sign(log2FoldChange)) %>%
  select(gene, sign, log2FoldChange, padj, baseMean)

cat(sprintf("G_tox: %d genes (padj < 0.05, |LFC| >= 1, baseMean >= 10)\n", nrow(g_tox)))
cat(sprintf("  %d up, %d down in toxic vs PBS\n",
            sum(g_tox$sign > 0), sum(g_tox$sign < 0)))

# ==============================================================================
# Step 3b: GSEA on pooled toxic vs PBS (Hallmark, GO:BP, Reactome)
# ==============================================================================

cat("\n=== GSEA (fgsea, MSigDB Hallmark + GO:BP + Reactome) ===\n")

# Rank metric: signed -log10(p), matching the prior toxseq GSEA convention.
ranks_df <- res_df %>%
  filter(!is.na(pvalue), !is.na(gene), gene != "", baseMean >= 10) %>%
  group_by(gene) %>%
  slice_min(pvalue, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(rank_metric = -log10(pvalue) * sign(log2FoldChange)) %>%
  arrange(desc(rank_metric))

gene_ranks <- setNames(ranks_df$rank_metric, ranks_df$gene)
cat(sprintf("Ranked %d genes (range %.2f to %.2f)\n",
            length(gene_ranks), min(gene_ranks), max(gene_ranks)))

# Human MSigDB mapped to mouse orthologs (matches prior toxseq GSEA).
msigdb_sets <- msigdbr(species = "Mus musculus")

build_pathways <- function(coll, subcoll = NULL) {
  x <- msigdb_sets %>% filter(gs_collection == coll)
  if (!is.null(subcoll)) x <- x %>% filter(gs_subcollection == subcoll)
  split(x$gene_symbol, x$gs_name)
}

hallmark_sets <- build_pathways("H")
gobp_sets     <- build_pathways("C5", "GO:BP")
reactome_sets <- build_pathways("C2", "CP:REACTOME")

cat(sprintf("  Hallmark: %d sets | GO:BP: %d | Reactome: %d\n",
            length(hallmark_sets), length(gobp_sets), length(reactome_sets)))

run_gsea <- function(pathways, label) {
  set.seed(42)
  res <- fgsea(pathways = pathways, stats = gene_ranks,
               minSize = 15, maxSize = 500, nPermSimple = 10000) %>%
    as_tibble() %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) %>%
    arrange(padj, pval)
  cat(sprintf("  %s: %d sig at padj < 0.05 (%d up, %d down)\n", label,
              sum(res$padj < 0.05, na.rm = TRUE),
              sum(res$padj < 0.05 & res$NES > 0, na.rm = TRUE),
              sum(res$padj < 0.05 & res$NES < 0, na.rm = TRUE)))
  res
}

gsea_hallmark <- run_gsea(hallmark_sets, "Hallmark")
gsea_gobp     <- run_gsea(gobp_sets,     "GO:BP   ")
gsea_reactome <- run_gsea(reactome_sets, "Reactome")

dir.create("transcriptomics/24-0432/output/gsea", showWarnings = FALSE, recursive = TRUE)

# Full results
write_tsv(gsea_hallmark, "transcriptomics/24-0432/output/gsea/gsea_24_toxic_vs_pbs_hallmark.tsv")
write_tsv(gsea_gobp,     "transcriptomics/24-0432/output/gsea/gsea_24_toxic_vs_pbs_gobp.tsv.gz")
write_tsv(gsea_reactome, "transcriptomics/24-0432/output/gsea/gsea_24_toxic_vs_pbs_reactome.tsv.gz")

# Filtered convenience tables (matches prior layout: hallmark_sig, gobp_top100, reactome_top100)
write_tsv(gsea_hallmark %>% filter(padj < 0.05),
          "transcriptomics/24-0432/output/gsea/gsea_24_toxic_vs_pbs_hallmark_sig.tsv")
write_tsv(gsea_gobp %>% filter(padj < 0.05) %>% head(100),
          "transcriptomics/24-0432/output/gsea/gsea_24_toxic_vs_pbs_gobp_top100.tsv")
write_tsv(gsea_reactome %>% filter(padj < 0.05) %>% head(100),
          "transcriptomics/24-0432/output/gsea/gsea_24_toxic_vs_pbs_reactome_top100.tsv")

cat("\nTop 10 Hallmark pathways:\n")
print(gsea_hallmark %>% select(pathway, padj, NES, size) %>% head(10))

# ==============================================================================
# Step 3c: Paraspeckle enrichment (custom gene set)
# ==============================================================================
# Paraspeckle-associated gene set curated per Fox et al., Trends Biochem Sci
# 43:124-135 (2018), https://www.cell.com/trends/biochemical-sciences/fulltext/S0968-0004(17)30224-4
# Mouse symbols are standard-capitalization orthologs of the published human
# list; BRG1 maps to Smarca4. Neat1 and Malat1 are lncRNAs included manually
# (not recovered by protein-ortholog mapping but central to paraspeckle biology).

paraspeckle_genes <- c(
  # Core structural components
  "Neat1", "Nono", "Sfpq", "Pspc1", "Rbm14",
  # RNA-binding proteins enriched in paraspeckles
  "Fus", "Tardbp", "Matr3", "Hnrnpk", "Hnrnpm", "Hnrnph1",
  "Hnrnph3", "Hnrnpa1", "Cpsf6", "Dazap1", "Smarca4",
  # Additional paraspeckle-associated factors
  "Ewsr1", "Taf15", "Rbmx", "Srsf10", "Celf6",
  "Hnrnpf", "Hnrnpul1", "Hnrnpa2b1", "Nudt21",
  "Malat1"
)

cat(sprintf("\n=== Paraspeckle enrichment (%d curated genes) ===\n",
            length(paraspeckle_genes)))

# ---- fgsea with the custom set against the existing ranked list ----
set.seed(42)
gsea_paraspeckle <- fgsea(
  pathways    = list(PARASPECKLE_COMPONENTS = paraspeckle_genes),
  stats       = gene_ranks,
  minSize     = 5,
  maxSize     = 500,
  nPermSimple = 10000
) %>%
  as_tibble() %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))

cat("fgsea result:\n")
print(gsea_paraspeckle %>% select(pathway, size, NES, pval, padj, leadingEdge))

# ---- Fisher's exact: over-representation among DE genes ----
# Universe = all genes retained by DESeq2 (res_df already drops NA padj).
universe_genes <- res_df$gene
sig_up   <- res_df %>% filter(padj < 0.05, log2FoldChange > 0) %>% pull(gene)
sig_down <- res_df %>% filter(padj < 0.05, log2FoldChange < 0) %>% pull(gene)
sig_all  <- res_df %>% filter(padj < 0.05) %>% pull(gene)

fisher_para <- function(query_genes, query_name) {
  a <- length(intersect(query_genes, paraspeckle_genes))
  b <- length(setdiff(paraspeckle_genes, query_genes))
  c <- length(setdiff(query_genes, paraspeckle_genes))
  d <- length(setdiff(universe_genes, union(query_genes, paraspeckle_genes)))
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
  tibble(
    query             = query_name,
    query_size        = length(query_genes),
    set_size_universe = length(intersect(paraspeckle_genes, universe_genes)),
    overlap           = a,
    odds_ratio        = unname(ft$estimate),
    p_value           = ft$p.value
  )
}

fisher_results <- bind_rows(
  fisher_para(sig_all,  "sig_all"),
  fisher_para(sig_up,   "sig_up"),
  fisher_para(sig_down, "sig_down")
) %>%
  mutate(padj = p.adjust(p_value, method = "BH"))

cat("Fisher's exact (one-sided, over-representation):\n")
print(fisher_results)

# ---- Per-gene DE table for the paraspeckle set ----
paraspeckle_de <- res_df %>%
  filter(gene %in% paraspeckle_genes) %>%
  select(gene, baseMean, log2FoldChange, lfcSE, pvalue, padj) %>%
  arrange(padj)

cat(sprintf("\nParaspeckle genes with DESeq2 results (%d / %d curated):\n",
            nrow(paraspeckle_de), length(paraspeckle_genes)))
print(paraspeckle_de)

write_tsv(gsea_paraspeckle, "transcriptomics/24-0432/output/gsea/paraspeckle_fgsea.tsv")
write_tsv(fisher_results,   "transcriptomics/24-0432/output/gsea/paraspeckle_fisher.tsv")
write_tsv(paraspeckle_de,   "transcriptomics/24-0432/output/gsea/paraspeckle_de_genes.tsv")

# ==============================================================================
# Step 4: VST + UMAP on G_tox genes
# ==============================================================================

cat("\n=== VST normalization ===\n")

# VST and UMAP on all samples (including modified) to preserve embedding geometry
meta <- meta_24
meta <- meta[match(colnames(count_mat), meta$sample_id), ]
rownames(meta) <- meta$sample_id

dds_all <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = meta,
  design    = ~ 1
)
vst_all <- vst(dds_all, blind = TRUE)
vst_mat <- assay(vst_all)

meta$condition <- meta$aso_class
print(table(meta$condition, meta$inoculum))

cat("\n=== UMAP ===\n")

g_tox_in_vst <- intersect(g_tox$gene, rownames(vst_mat))
vst_umap <- vst_mat[g_tox_in_vst, ]

set.seed(42)
umap_result <- umap(t(vst_umap), n_neighbors = 30, min_dist = 0.8, spread = 2, n_components = 2)
meta$umap1 <- umap_result[, 1]
meta$umap2 <- umap_result[, 2]

cat(sprintf("UMAP computed on %d G_tox genes\n", length(g_tox_in_vst)))

# ==============================================================================
# Step 5: G_tox continuous score
# ==============================================================================

cat("\n=== Computing G_tox scores ===\n")

# ---- Load pre-computed per-ASO contrasts ----
contrasts_24_raw <- read_tsv("analytic/all_contrasts_0432.tsv.gz",
                             show_col_types = FALSE)
# Keep ASO-vs-PBS contrasts for both inoculums (exclude PBS-vs-PBS)
contrasts_24 <- contrasts_24_raw %>%
  filter(grepl("_PBS$", group2), !grepl("_PBS$", group1)) %>%
  mutate(contrast_id = paste0(group1, "_vs_", group2))

# ---- G_tox continuous L2FC score ----
g_tox_sign <- g_tox %>% select(gene, sign)

compute_gtox_continuous <- function(df) {
  df %>%
    inner_join(g_tox_sign, by = "gene") %>%
    group_by(contrast_id) %>%
    summarise(
      gtox_score = sum(log2fold_change * sign, na.rm = TRUE),
      n_gtox = n(),
      .groups = "drop"
    )
}

gtox_cont_24 <- compute_gtox_continuous(contrasts_24)

# ---- Merge, annotate, exclude modified ----
scores_24 <- gtox_cont_24 %>%
  mutate(
    parts    = str_match(contrast_id, "^(RML|CBH)_(.+)_vs_"),
    inoculum = parts[, 2],
    ionis_id = parts[, 3]
  ) %>%
  select(-parts) %>%
  left_join(aso_lookup_24 %>% select(ionis_id, aso_name), by = "ionis_id") %>%
  mutate(
    condition = case_when(
      ionis_id %in% TOXIC_ASOS     ~ "toxic",
      ionis_id == ACTIVE_ASO2      ~ "toxic",
      ionis_id %in% MODIFIED_ASOS  ~ "modified",
      ionis_id %in% TOLERATED_ASOS ~ "tolerated",
      TRUE                         ~ "other"
    ),
    label = paste0(aso_name, " (", inoculum, ")")
  ) %>%
  filter(condition != "modified") %>%
  arrange(gtox_score)

cat("\n--- Study 24 G_tox Scores ---\n")
scores_24 %>%
  select(label, inoculum, gtox_score, condition) %>%
  print(n = Inf)

# ==============================================================================
# Step 6: Figure — UMAP + G_tox score dot plot
# ==============================================================================

cat("\n=== Generating figure ===\n")

condition_colors <- c(
  toxic     = "tomato",
  tolerated = "cornflowerblue",
  PBS       = "gray60"
)

png("transcriptomics/24-0432/output/fig1.png", width = 14, height = 6, units = "in", res = 300)
layout(matrix(1:2, nrow = 1), widths = c(1, 1))

# -- Panel A: UMAP --
umap_df <- meta %>% filter(condition != "modified")
umap_pad <- 1.5
umap_xlim <- c(floor(min(umap_df$umap1)) - umap_pad,
               ceiling(max(umap_df$umap1)) + umap_pad)
umap_ylim <- c(floor(min(umap_df$umap2)) - umap_pad,
               ceiling(max(umap_df$umap2)) + umap_pad)

par(mar = c(4, 9, 3, 1))
plot(NA, xlim = umap_xlim, ylim = umap_ylim,
     xlab = "", ylab = "", axes = FALSE, frame.plot = FALSE)

for (cond in unique(umap_df$condition)) {
  for (inoc in unique(umap_df$inoculum)) {
    sub <- umap_df %>% filter(condition == cond, inoculum == inoc)
    if (nrow(sub) == 0) next
    pch_val <- if (inoc == "RML") 16 else 1
    points(sub$umap1, sub$umap2,
           col = condition_colors[cond], pch = pch_val, cex = 1.2,
           lwd = if (inoc == "CBH") 2 else 1)
  }
}

umap_xticks <- pretty(umap_xlim)
umap_yticks <- pretty(umap_ylim)
axis(1, at = umap_xticks, lwd = 0.5)
axis(2, at = umap_yticks, las = 1, lwd = 0.5)
mtext("UMAP 1", side = 1, line = 2.5)
mtext("UMAP 2", side = 2, line = 2.5)
mtext("A", side = 3, adj = 0, font = 2, cex = 1.2, line = 1)

legend("topright",
       legend = c("toxic", "tolerated", "PBS", "RML (solid)", "CBH (hollow)"),
       col = c(condition_colors[c("toxic", "tolerated", "PBS")], "black", "black"),
       pch = c(16, 16, 16, 16, 1),
       pt.lwd = c(1, 1, 1, 1, 2),
       bty = "n", cex = 0.85)

# -- Panel B: G_tox Score --
n_asos <- nrow(scores_24)
score_range <- range(scores_24$gtox_score)
score_xticks <- pretty(score_range, n = 5)
score_xlim <- range(score_xticks)

par(mar = c(4, 9, 3, 1))
plot(NA, xlim = score_xlim, ylim = c(0.5, n_asos + 0.5),
     xlab = "", ylab = "", axes = FALSE, frame.plot = FALSE)

for (i in seq_len(n_asos)) {
  row <- scores_24[i, ]
  pch_val <- if (row$inoculum == "RML") 16 else 1
  points(row$gtox_score, i,
         col = condition_colors[row$condition], pch = pch_val, cex = 1.4,
         lwd = if (row$inoculum == "CBH") 2 else 1)
}

axis(1, at = score_xticks, lwd = 0.5)
axis(2, at = seq_len(n_asos), labels = scores_24$label,
     las = 2, cex.axis = 0.85, lwd = 0.5, tick = FALSE)
mtext(expression(Sigma~log2FC~G[tox]), side = 1, line = 2.5)
mtext("B", side = 3, adj = 0, font = 2, cex = 1.2, line = 1)
abline(v = 0, lty = 3, col = "grey50")

dev.off()
cat("Saved: transcriptomics/24-0432/output/fig1.png\n")

# ---- Fig 2: Volcano plot ----
astro_genes <- read_csv("analytic/astrocyte_genes.csv", show_col_types = FALSE) %>%
  clean_names() %>%
  filter(astrocyte_average >= 0.1,
         astrocyte_log2_fold_change >= 2,
         astrocyte_p_value < 0.05) %>%
  pull(feature_name)

volcano_df <- res_df %>%
  mutate(
    is_astro = gene %in% astro_genes,
    color_cat = ifelse(is_astro, "Astrocyte", "Other"),
    neg_log10p = -log10(pvalue)
  )

astro_label <- volcano_df %>%
  filter(is_astro, padj < 0.05) %>%
  slice_min(padj, n = 20)

p_volcano <- ggplot(volcano_df,
                    aes(x = log2FoldChange, y = neg_log10p, color = color_cat)) +
  geom_point(alpha = 0.5, size = 1) +
  ggrepel::geom_text_repel(
    data = astro_label, aes(label = gene),
    size = 3.5, show.legend = FALSE, max.overlaps = 20,
    min.segment.length = 0.2, segment.color = "grey50",
    force = 2, force_pull = 0.5
  ) +
  scale_color_manual(values = c(
    "Astrocyte" = "#CC6576",
    "Other"     = "grey70"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey40") +
  labs(
    x     = expression(log[2]~Fold~Change~"(Toxic / PBS)"),
    y     = expression(-log[10]~italic(p)),
    color = NULL,
    title = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = c(0.85, 0.85),
        legend.background = element_rect(fill = "white", color = NA))

ggsave("transcriptomics/24-0432/output/fig2.png", p_volcano, width = 8, height = 6, dpi = 300)
cat("Saved: transcriptomics/24-0432/output/fig2.png\n")

# ---- Export tables ----
write_csv(scores_24 %>%
            select(contrast_id, aso = ionis_id, aso_name, condition, inoculum, label, gtox_score),
          "transcriptomics/24-0432/output/dge_scores_24.csv")

write_csv(g_tox %>%
            select(gene, sign, log2FoldChange, padj, baseMean),
          "transcriptomics/24-0432/output/g_tox_genes.csv")

cat("Saved: transcriptomics/24-0432/output/dge_scores_24.csv\n")
cat("Saved: transcriptomics/24-0432/output/g_tox_genes.csv\n")

cat("\n=== Done ===\n")
