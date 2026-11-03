# =============================================================================
# Taxonomic attribution of DE KOs and pathway-level responses (GSEA)
#
#   A-B  family-level taxonomic attribution of the significant KOs, per contrast,
#        for Bacteria and Fungi
#   C-D  GSEA normalised enrichment scores of selected KEGG pathways
#
# GSEA uses the official KEGG KO -> map association (ko00001); the data file
# gsea_nes_selected.tsv is built by prepare_gsea_nes_heatmap.py from
# results/gsea/gsea_official_all.tsv. Global maps (e.g. Biosynthesis of amino
# acids ko01230, Carbon metabolism ko01200) have no KO list in ko00001 and are
# not tested.
#
# Run from this script's folder:  Rscript plot_taxonomic_attribution_and_gsea.R
# (FIG_DATA = folder holding de_family_attribution_plot.tsv and gsea_nes_selected.tsv)
# (FIG_DIR may point at another folder holding _style.R and data/)
# =============================================================================
HERE <- Sys.getenv("FIG_DIR", ".")
source(file.path(HERE, "_style.R"))

D   <- Sys.getenv("FIG_DATA", file.path(HERE, "data"))
OUT <- HERE

E_COL_MM <- 12.0; E_ROW_MM <- 5.6

# Subtitles 10 pt and panel tags 12 pt
add_tag <- function(g, lab) {
  g2 <- gtable_add_rows(g, unit(5, "mm"), pos = 0)
  gtable_add_grob(g2,
    textGrob(lab, x = unit(0.4, "mm"), y = unit(0.2, "mm"), hjust = 0, vjust = 0,
             gp = gpar(fontsize = 12, fontface = "bold", fontfamily = FONT)),
    t = 1, l = 1, r = ncol(g2), name = "tag", clip = "off")
}

# =============================================================================
# A - taxonomic attribution (family level)
# =============================================================================
DD <- read.delim(file.path(D, "de_family_attribution_plot.tsv")) |>
  mutate(contrast_label = factor(CONTRAST_LEVELS[match(contrast, CONTRAST_KEYS)],
                                 levels = CONTRAST_LEVELS))

fam_levels <- function(dom) {
  f <- DD |> filter(domain == dom, family != "No DE KOs") |>
    group_by(family) |> summarise(n = sum(n), .groups = "drop") |>
    arrange(desc(n)) |> pull(family)
  c(setdiff(f, "Other families"), intersect("Other families", f))
}
# alternate hue by stacking rank so neighbouring bands never share a hue
fam_cols <- function(levs) {
  named <- setdiff(levs, "Other families")
  pal <- c("#6B3FA0", "#1B7837", "#C2A5CF", "#A6DBA0", "#40004B", "#00441B",
           "#9970AB", "#5AAE61", "#E7D4E8", "#D9F0D3", "#762A83", "#238B45")
  setNames(c(rep(pal, length.out = max(length(named), 1))[seq_along(named)], "grey82"),
           c(named, "Other families"))
}
panelA <- function(dom, leg_ncol) {
  d <- filter(DD, domain == dom, family != "No DE KOs")
  levs <- fam_levels(dom); d$family <- factor(d$family, levels = levs)
  empty <- DD |> filter(domain == dom, family == "No DE KOs")
  nlab  <- d |> group_by(contrast_label) |>
    summarise(n_total = first(n_total), .groups = "drop")
  ggplot(d, aes(contrast_label, pct, fill = family)) +
    geom_col(width = 0.66, colour = NA) +
    # bare count: "n=445" overruns the ~6.3 mm the bar pitch allows
    geom_text(data = nlab, aes(contrast_label, 101, label = n_total),
              inherit.aes = FALSE, vjust = 0, size = BASE / .pt, family = FONT,
              colour = "grey25") +
    geom_text(data = empty, aes(contrast_label, 50, label = "no DE"),
              inherit.aes = FALSE, angle = 90, size = BASE / .pt,
              family = FONT, colour = "grey45") +
    geom_vline(xintercept = c(3.5, 6.5), linetype = "dashed",
               colour = "grey45", linewidth = 0.25) +
    scale_fill_manual(values = fam_cols(levs), name = NULL,
                      guide = guide_legend(ncol = leg_ncol, byrow = FALSE)) +
    scale_x_discrete(limits = CONTRAST_LEVELS, drop = FALSE) +
    scale_y_continuous(limits = c(0, 109), breaks = seq(0, 100, 25),
                       expand = expansion(mult = c(0, 0))) +
    labs(subtitle = dom, x = NULL, y = "Significant KOs\nassigned to family (%)") +
    theme_fig() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = BASE),
          legend.position = "bottom", legend.key.size = unit(0.24, "cm"),
          legend.box.margin = margin(0, 0, 0, 0),
          plot.subtitle = element_text(size = 10, colour = if (dom == "Bacteria") BLUE else PURPLE))
}
gA <- cbind(add_tag(set_panel_mm(panelA("Bacteria", 2), 57, 34), "A"),
            add_tag(set_panel_mm(panelA("Fungi", 1),    57, 34), "B"), size = "max")

# =============================================================================
# B - GSEA v2 (deduplicated features, BH over the full pathway universe)
# =============================================================================
E <- read.delim(file.path(D, "gsea_nes_selected.tsv"), check.names = FALSE,
                fileEncoding = "UTF-8")
COL_LEVELS <- c("G24/C24", "G48/C48", "G72/C72",
                "\u0394 effect\n72\u221224 h")     # \u escapes: locale-proof
E$col <- factor(E$col, levels = unique(E$col))
levels(E$col) <- COL_LEVELS[seq_along(levels(E$col))]

E_LIM <- 2.3
panelB <- function(dom, row_order, relabel = NULL) {
  d <- filter(E, domain == dom)
  d$lab_y <- as.character(d$pathway_name)
  if (!is.null(relabel)) {
    hit <- d$lab_y %in% names(relabel); d$lab_y[hit] <- relabel[d$lab_y[hit]]
    row_order <- ifelse(row_order %in% names(relabel), relabel[row_order], row_order)
  }
  d$lab_y <- factor(d$lab_y, levels = rev(row_order))
  # the text colour goes in the DATA and is mapped with an identity scale;
  # a bare vector passed to colour= is not guaranteed to stay in row order
  d$ink <- on_viridis(d$NES, -E_LIM, E_LIM)
  ggplot(d, aes(col, lab_y, fill = NES)) +
    geom_tile() +
    geom_text(aes(label = paste0(sprintf("%.1f", NES), mark), colour = ink),
              size = BASE / .pt, family = FONT, fontface = "bold") +
    scale_colour_identity() +
    scale_fill_gradientn(colours = VIR, limits = c(-E_LIM, E_LIM),
                         oob = scales::squish, name = "NES",
                         breaks = c(-2, -1, 0, 1, 2)) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(subtitle = dom, x = NULL, y = NULL) +
    theme_heat() +
    theme(axis.text.y = element_text(size = BASE, colour = "black"),
          axis.text.x.top = element_text(size = BASE, face = "bold"),
          plot.subtitle = element_text(hjust = 0, size = 10,
                            colour = if (dom == "Bacteria") BLUE else PURPLE),
          legend.position = "none")
}
ordBb <- c("Glycolysis / Gluconeogenesis", "Pentose phosphate pathway",
           "Nitrogen metabolism", "Pyruvate metabolism", "Citrate cycle (TCA cycle)", "Bacterial chemotaxis",
           "Flagellar assembly", "Ribosome", "Peptidoglycan biosynthesis")
ordBf <- c("Ribosome", "Glycolysis / Gluconeogenesis", "Pyruvate metabolism", "Citrate cycle (TCA cycle)",
           "Valine, leucine and isoleucine degradation", "Butanoate metabolism",
           "beta-Alanine metabolism")
RELAB <- c("Valine, leucine and isoleucine degradation" = "Val/Leu/Ile degradation",
           "beta-Alanine metabolism" = "\u03b2-Alanine metabolism")
# the NES key rides on the fungal sub-panel so the two heatmaps can simply be
# column-bound; a free-floating legend grob has nothing to align to
stopifnot(setequal(ordBb, unique(E$pathway_name[E$domain == "Bacteria"])),
          setequal(ordBf, unique(E$pathway_name[E$domain == "Fungi"])))
gB <- cbind(
  add_tag(set_panel_mm(panelB("Bacteria", ordBb), 4 * E_COL_MM, length(ordBb) * E_ROW_MM), "C"),
  add_tag(set_panel_mm(panelB("Fungi", ordBf, RELAB) +
                 theme(legend.position = "right",
                       legend.key.width  = unit(0.20, "cm"),
                       legend.key.height = unit(0.85, "cm")),
               4 * E_COL_MM, length(ordBf) * E_ROW_MM), "D"),
  size = "max")

# =============================================================================
# assemble
# =============================================================================
PAD <- 2
W <- max(gw(gA), gw(gB))

fig <- wrap_elements(full = pad_right(gA, W)) /
       wrap_elements(full = pad_right(gB, W)) +
  plot_layout(heights = unit(c(gh(gA), gh(gB)) + PAD, "mm"),
              widths  = unit(W + PAD, "mm")) +
  plot_annotation(theme = theme(plot.margin = margin(2, 2, 2, 2))) &
  theme(plot.margin = margin(0, 0, 0, 0))

W_MM <- W + PAD + 4
H_MM <- gh(gA) + gh(gB) + 2 * PAD + 4
save_fig(fig, file.path(OUT, "taxonomic_attribution_and_gsea"), W_MM, H_MM)

cat(sprintf("  A %.0f x %.0f | B %.0f x %.0f mm\n",
            gw(gA), gh(gA), gw(gB), gh(gB)))
cat("panel B rows:", length(ordBb), "bacterial +", length(ordBf), "fungal\n")
cat("smallest font on the figure:", BASE, "pt\n")
invisible(dev.off())
if (file.exists(file.path(OUT, "Rplots.pdf"))) unlink(file.path(OUT, "Rplots.pdf"))
