# =============================================================================
# KO-level differential expression -- glucose incubation metatranscriptome
# Implements DE_analysis_scope.pdf (2026-08-17) exactly.
#
# Bacteria and Fungi are analysed as two INDEPENDENT experiments: separate count
# matrices, separate size factors, separate dispersion estimation, separate BH
# correction. Archaea are excluded (Scope B). No KO is filtered on TRI, prior
# significance or taxonomic family.
#
# Design: ~ treatment * timepoint   (treatment: control, glucose;
#                                    timepoint: 24h, 48h, 72h; n = 3 per cell)
#
# 12 contrasts per domain, all extracted from that ONE fitted model as numeric
# contrast vectors over its six coefficients -- never by refitting on subsets,
# so dispersion is estimated once from all 18 samples:
#
#   GLUCOSE EFFECT      G24 vs C24, G48 vs C48, G72 vs C72
#   TEMPORAL (control)  C48 vs C24, C72 vs C24, C72 vs C48
#   TEMPORAL (glucose)  G48 vs G24, G72 vs G24, G72 vs G48
#   INTERACTION         (G48-C48)-(G24-C24), (G72-C72)-(G24-C24),
#                       (G72-C72)-(G48-C48)
#
# BH is applied by DESeq2 WITHIN each contrast, within each domain. Bacteria and
# Fungi are never pooled into one correction, and the 12 contrasts are never
# pooled into a global one -- both per the scope.
#
# Inputs (small, prebuilt by 02_count_matrices/build_ko_domain_matrices.py -- the 108 MiB
# counts_ko_taxonomy_full.tsv is NOT read here):
#   counts_ko_bacteria.tsv, counts_ko_fungi.tsv        KO x 18 raw counts
#   ko_top_family_bacteria.tsv, ko_top_family_fungi.tsv
#   ko_to_pathway.tsv, sample_metadata.tsv
#
# Run: Rscript deseq2_ko_by_domain.R
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
})

PROJECT_DIR <- Sys.getenv("MTX_ROOT", ".")
IN_DIR  <- file.path(PROJECT_DIR, "results", "ko_matrices")
OUT_DIR <- file.path(PROJECT_DIR, "results", "deseq2")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

LFC_HIGHLIGHT <- 1      # |log2FC| >= 1 flags a biologically larger response.
ALPHA         <- 0.05   # padj threshold -- the statistical criterion.

# --- design -------------------------------------------------------------------
meta <- read.delim(file.path(IN_DIR, "sample_metadata.tsv"), stringsAsFactors = FALSE)
rownames(meta) <- meta$sample_id
meta$treatment <- factor(meta$treatment, levels = c("control", "glucose"))
meta$timepoint <- factor(meta$timepoint, levels = c("24h", "48h", "72h"))
stopifnot(nrow(meta) == 18, all(table(meta$treatment, meta$timepoint) == 3))

# Coefficient order of ~ treatment * timepoint with those reference levels:
#   1 Intercept
#   2 treatment_glucose_vs_control
#   3 timepoint_48h_vs_24h
#   4 timepoint_72h_vs_24h
#   5 treatmentglucose.timepoint48h
#   6 treatmentglucose.timepoint72h
COEFS <- c("Intercept", "treatment_glucose_vs_control",
           "timepoint_48h_vs_24h", "timepoint_72h_vs_24h",
           "treatmentglucose.timepoint48h", "treatmentglucose.timepoint72h")

# Each vector is the linear combination of those coefficients that equals the
# difference in the comment. Verified against the model matrix by
# check_contrasts.R before this script was ever run.
CONTRASTS <- list(
  # ---- glucose effect at each timepoint ----
  list(id = "glucose_G24_vs_C24", family = "glucose effect", label = "G24 vs C24",
       v = c(0, 1, 0, 0, 0, 0)),
  list(id = "glucose_G48_vs_C48", family = "glucose effect", label = "G48 vs C48",
       v = c(0, 1, 0, 0, 1, 0)),
  list(id = "glucose_G72_vs_C72", family = "glucose effect", label = "G72 vs C72",
       v = c(0, 1, 0, 0, 0, 1)),
  # ---- temporal trajectory within control ----
  list(id = "time_C48_vs_C24", family = "temporal (control)", label = "C48 vs C24",
       v = c(0, 0, 1, 0, 0, 0)),
  list(id = "time_C72_vs_C24", family = "temporal (control)", label = "C72 vs C24",
       v = c(0, 0, 0, 1, 0, 0)),
  list(id = "time_C72_vs_C48", family = "temporal (control)", label = "C72 vs C48",
       v = c(0, 0, -1, 1, 0, 0)),
  # ---- temporal trajectory within glucose ----
  list(id = "time_G48_vs_G24", family = "temporal (glucose)", label = "G48 vs G24",
       v = c(0, 0, 1, 0, 1, 0)),
  list(id = "time_G72_vs_G24", family = "temporal (glucose)", label = "G72 vs G24",
       v = c(0, 0, 0, 1, 0, 1)),
  list(id = "time_G72_vs_G48", family = "temporal (glucose)", label = "G72 vs G48",
       v = c(0, 0, -1, 1, -1, 1)),
  # ---- interaction: does the glucose effect change with time? ----
  list(id = "interaction_48_vs_24", family = "interaction",
       label = "(G48-C48) - (G24-C24)", v = c(0, 0, 0, 0, 1, 0)),
  list(id = "interaction_72_vs_24", family = "interaction",
       label = "(G72-C72) - (G24-C24)", v = c(0, 0, 0, 0, 0, 1)),
  list(id = "interaction_72_vs_48", family = "interaction",
       label = "(G72-C72) - (G48-C48)", v = c(0, 0, 0, 0, -1, 1))
)

# --- annotation ---------------------------------------------------------------
pathway <- read.delim(file.path(IN_DIR, "ko_to_pathway.tsv"), stringsAsFactors = FALSE)

# KO -> gene name. The project's ko_names_lookup.tsv covers only ~60 KOs, so the
# full dictionary is pulled from the KEGG REST list endpoint and cached beside
# the inputs. If there is no network the run continues with names left blank --
# the DE numbers do not depend on it.
ko_names_file <- file.path(IN_DIR, "ko_names_full.tsv")
if (file.exists(ko_names_file)) {
  ko_names <- read.delim(ko_names_file, stringsAsFactors = FALSE)
  message(sprintf("KO names: %d from cache", nrow(ko_names)))
} else {
  ko_names <- tryCatch({
    x <- read.delim("https://rest.kegg.jp/list/ko", header = FALSE,
                    quote = "", stringsAsFactors = FALSE,
                    col.names = c("KO", "gene_name"))
    x$KO <- sub("^ko:", "", x$KO)
    write.table(x, ko_names_file, sep = "\t", row.names = FALSE, quote = FALSE)
    message(sprintf("KO names: %d fetched from KEGG and cached", nrow(x)))
    x
  }, error = function(e) {
    message("KO names: KEGG unreachable (", conditionMessage(e),
            ") -- falling back to the local lookup")
    f <- file.path(IN_DIR, "ko_names_lookup.tsv")
    if (file.exists(f)) read.delim(f, stringsAsFactors = FALSE)
    else data.frame(KO = character(0), gene_name = character(0))
  })
}

# --- per domain ---------------------------------------------------------------
all_res <- list()
summary_rows <- list()
dds_list <- list()

for (dom in c("Bacteria", "Fungi")) {
  tag <- tolower(dom)
  message(sprintf("\n=================== %s ===================", dom))

  cts <- read.delim(file.path(IN_DIR, sprintf("counts_ko_%s.tsv", tag)),
                    row.names = 1, check.names = FALSE)
  cts <- as.matrix(cts[, rownames(meta), drop = FALSE])
  storage.mode(cts) <- "integer"
  message(sprintf("  %d KOs x %d samples | %s reads",
                  nrow(cts), ncol(cts), format(sum(cts), big.mark = ",")))

  topfam <- read.delim(file.path(IN_DIR, sprintf("ko_top_family_%s.tsv", tag)),
                       stringsAsFactors = FALSE)

  dds <- DESeqDataSetFromMatrix(cts, colData = meta, design = ~ treatment * timepoint)
  dds <- DESeq(dds, quiet = TRUE)
  dds_list[[dom]] <- dds

  # The contrast vectors above are positional; assert the fitted model really has
  # the coefficient order they assume, so a DESeq2 change can never silently
  # remap them onto the wrong terms.
  stopifnot(identical(resultsNames(dds),
                      c("Intercept", "treatment_glucose_vs_control",
                        "timepoint_48h_vs_24h", "timepoint_72h_vs_24h",
                        "treatmentglucose.timepoint48h",
                        "treatmentglucose.timepoint72h")))

  sf <- sizeFactors(dds)
  message(sprintf("  size factors %.3f .. %.3f", min(sf), max(sf)))

  for (ct in CONTRASTS) {
    r <- results(dds, contrast = ct$v, alpha = ALPHA)
    df <- as.data.frame(r)
    df$KO <- rownames(df)
    df$domain <- dom
    df$contrast <- ct$id
    df$contrast_label <- ct$label
    df$contrast_family <- ct$family

    df <- merge(df, ko_names, by = "KO", all.x = TRUE)
    df <- merge(df, pathway, by = "KO", all.x = TRUE)
    df <- merge(df, topfam[, c("KO", "top_family", "top_family_pct")],
                by = "KO", all.x = TRUE)

    df$significant <- !is.na(df$padj) & df$padj < ALPHA
    df$significant_lfc1 <- df$significant & abs(df$log2FoldChange) >= LFC_HIGHLIGHT

    df <- df[order(df$padj, df$pvalue), c(
      "KO", "gene_name", "domain", "contrast", "contrast_label", "contrast_family",
      "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj",
      "significant", "significant_lfc1", "top_family", "top_family_pct",
      "KEGG_Pathway")]

    write.table(df, file.path(OUT_DIR, sprintf("de_%s_%s.tsv", tag, ct$id)),
                sep = "\t", row.names = FALSE, quote = FALSE)
    all_res[[paste(dom, ct$id)]] <- df

    n_sig <- sum(df$significant)
    summary_rows[[paste(dom, ct$id)]] <- data.frame(
      domain = dom, contrast = ct$id, contrast_label = ct$label,
      contrast_family = ct$family,
      n_tested = sum(!is.na(df$padj)), n_ko_in_matrix = nrow(df),
      n_significant = n_sig,
      n_up = sum(df$significant & df$log2FoldChange > 0),
      n_down = sum(df$significant & df$log2FoldChange < 0),
      n_significant_lfc1 = sum(df$significant_lfc1),
      min_padj = suppressWarnings(min(df$padj, na.rm = TRUE)))
    message(sprintf("  %-22s %-24s %5d significant (%d up / %d down)",
                    ct$id, ct$label, n_sig,
                    sum(df$significant & df$log2FoldChange > 0),
                    sum(df$significant & df$log2FoldChange < 0)))
  }
}

combined <- do.call(rbind, all_res)
write.table(combined, file.path(OUT_DIR, "de_all_contrasts_long.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

summ <- do.call(rbind, summary_rows)
rownames(summ) <- NULL
write.table(summ, file.path(OUT_DIR, "summary_de_by_contrast.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Normalised counts, for the DE heatmaps -- so the plots never have to
# re-derive them and can never disagree with the DE table.
for (dom in names(dds_list)) {
  nc <- counts(dds_list[[dom]], normalized = TRUE)
  write.table(data.frame(KO = rownames(nc), nc, check.names = FALSE),
              file.path(OUT_DIR, sprintf("normalized_counts_%s.tsv", tolower(dom))),
              sep = "\t", row.names = FALSE, quote = FALSE)
  saveRDS(dds_list[[dom]], file.path(OUT_DIR, sprintf("dds_%s.rds", tolower(dom))))
}

writeLines(
  c(sprintf("run finished: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("DESeq2 version: %s", as.character(utils::packageVersion("DESeq2"))),
    sprintf("design: ~ treatment * timepoint | alpha = %.2f | |log2FC| flag = %.1f",
            ALPHA, LFC_HIGHLIGHT),
    sprintf("filter applied upstream: total count >= 10 AND detected in >= 3 of 18 samples"),
    "", capture.output(sessionInfo())),
  file.path(OUT_DIR, "deseq2_session_info.txt"))

cat("\n==================== SUMMARY ====================\n")
print(summ[, c("domain", "contrast_label", "contrast_family",
               "n_significant", "n_up", "n_down", "n_significant_lfc1")],
      row.names = FALSE)
cat(sprintf("\nOutputs written to: %s\n", OUT_DIR))
