# =============================================================================
# Gene-level differential expression: DE counts per contrast and clustered heatmaps
#
#   A  number of DE KOs per contrast, Bacteria and Fungi
#   B  every bacterial KO significant in at least one contrast,
#      rows only clustered, no row labels, very thin rows
#   C  every fungal KO significant anywhere + functional description
#
# PANEL B IS DRAWN WITH geom_raster(), NOT geom_tile()
# With ~1000 rows a tile is ~0.13 mm tall. geom_tile emits rows x 6 separate
# rectangles, and their antialiased edges overlap into visible banding that is
# not in the data. geom_raster rasterises the matrix in one pass at device
# resolution, so what you see is the matrix. The cluster strip beside it stays
# geom_tile because it needs an explicit width.
#
# Panel C labels each fungal KO with the KEGG map on which ko00001 places it,
# or with its function when it is on no map (e.g. clpC, a Clp chaperone/protease).
#
# Run from this script's folder:  Rscript plot_de_counts_and_heatmaps.R
# (FIG_DATA = folder written by prepare_de_heatmap_clusters.py)
# =============================================================================
HERE <- Sys.getenv("FIG_DIR", ".")
source(file.path(HERE, "_style.R"))

D   <- Sys.getenv("FIG_DATA", file.path(HERE, "data"))
OUT <- HERE

B_COL_MM <- 4.0                     # bacterial heatmap column width
B_H_MM   <- 148                     # bacterial heatmap total height
C_COL_MM <- 6.0;  C_ROW_MM <- 5.6   # fungal heatmap tile pitch

# =============================================================================
# A - number of DE KOs per contrast
# =============================================================================
A <- read.delim(file.path(D, "de_counts_by_contrast.tsv")) |>
  mutate(contrast_label = factor(contrast_label, levels = CONTRAST_LEVELS),
         direction = factor(direction, levels = c("Upregulated", "Downregulated")))
totA <- A |> group_by(domain, contrast_label) |> summarise(n = sum(n), .groups = "drop")

# Horizontal bars. coord_flip() only swaps the rendering, so everything is
# still written in the pre-flip frame: contrast is the DISCRETE x, count the
# continuous y, separators stay geom_vline(). The discrete order is reversed
# because flipping puts the first level at the bottom.
panelA <- function(dom, cmax, brk) {
  d <- filter(A, domain == dom); t <- filter(totA, domain == dom)
  ggplot(d, aes(contrast_label, n, fill = direction)) +
    geom_col(width = 0.66, colour = NA) +
    geom_text(data = t, aes(contrast_label, n, label = n), inherit.aes = FALSE,
              hjust = -0.25, size = BASE / .pt, family = FONT, fontface = "bold") +
    geom_vline(xintercept = c(3.5, 6.5), linetype = "dashed",
               colour = "grey45", linewidth = 0.25) +
    scale_fill_manual(values = c(Upregulated = PURPLE, Downregulated = GREEN),
                      name = NULL) +
    scale_x_discrete(limits = rev(CONTRAST_LEVELS), drop = FALSE) +
    scale_y_continuous(limits = c(0, cmax), breaks = brk,
                       expand = expansion(mult = c(0, 0.11))) +
    coord_flip() +
    labs(subtitle = dom, y = "Number of DE KOs", x = NULL) +
    theme_fig() +
    theme(axis.text.y = element_text(size = BASE, hjust = 1),
          legend.position = "none",
          plot.subtitle = element_text(hjust = 0,
                            colour = if (dom == "Bacteria") BLUE else PURPLE))
}
gA <- rbind(
  set_panel_mm(panelA("Bacteria", 1000, seq(0, 1000, 250)), 104, 30),
  set_panel_mm(panelA("Fungi", 3.3, 0:3) +
                 theme(legend.position = "bottom",
                       legend.margin = margin(0, 0, 0, 0)), 104, 30),
  size = "max")

# =============================================================================
# B - all 1073 significant bacterial KOs
# =============================================================================
wB  <- read.delim(file.path(D, "de_heatmap_bacteria_z.tsv"), row.names = 1,
                  check.names = FALSE)[, GROUPS]
clB <- read.delim(file.path(D, "de_heatmap_bacteria_clusters.tsv"))
prof <- read.delim(file.path(D, "de_heatmap_cluster_profiles.tsv"))
nB <- nrow(wB)
stopifnot(identical(rownames(wB), clB$KO))     # order comes from the prep script

# row 1 of the file is the TOP of the heatmap, so y counts down from nB
lB <- as.data.frame(as.table(as.matrix(wB))); names(lB) <- c("KO", "group", "z")
lB$xi <- match(as.character(lB$group), GROUPS)
lB$yi <- nB - match(as.character(lB$KO), rownames(wB)) + 1
clB$yi <- nB - clB$row + 1

# Dendrogram (left) and cluster strip (right) both live in THIS plot, on the
# same continuous x axis as the matrix. Drawing either as a separate ggplot and
# aligning them afterwards is what makes rows drift; here they cannot.
#
# The tree is read from the prep script's export, not recomputed: it has to be
# the same tree the cluster assignments came from.
segB <- read.delim(file.path(D, "de_heatmap_dendrogram.tsv"))
hmax <- max(segB$height0, segB$height1)

DG_R <- -0.55; DG_L <- -3.65        # dendrogram occupies x in [DG_L, DG_R]
h2x  <- function(h) DG_R - (DG_R - DG_L) * (h / hmax)
segB$x0 <- h2x(segB$height0); segB$x1 <- h2x(segB$height1)
segB$y0 <- nB - segB$row0 + 1;  segB$y1 <- nB - segB$row1 + 1

CL_X  <- 7.05; CL_W <- 0.55         # cluster strip: the column after the matrix
CL_PAL <- c("#40004B", "#9970AB", "#C2A5CF", "#A6DBA0", "#1B7837",
            "#5AAE61", "#762A83", "#00441B")
ncl <- max(clB$cluster)
cl_cols <- setNames(CL_PAL[seq_len(ncl)], seq_len(ncl))
cl_mid <- clB |> group_by(cluster) |>
  summarise(yi = mean(yi), n = dplyr::n(), .groups = "drop")

XB_MIN <- -3.75; XB_MAX <- 8.5
pB <- ggplot() +
  geom_segment(data = segB, aes(x = x0, y = y0, xend = x1, yend = y1),
               linewidth = 0.12, colour = "grey35") +
  geom_raster(data = lB, aes(xi, yi, fill = z)) +
  geom_tile(data = clB, aes(CL_X, yi, fill = NULL),
            width = CL_W, height = 1,
            fill = cl_cols[as.character(clB$cluster)]) +
  geom_text(data = cl_mid, aes(CL_X + 0.45, yi, label = cluster),
            hjust = 0, size = BASE / .pt, family = FONT, fontface = "bold") +
  annotate("text", x = 1:6, y = nB * 1.012, label = GROUPS, angle = 90,
           hjust = 0, fontface = "bold", size = BASE / .pt, family = FONT) +
  annotate("text", x = CL_X, y = nB * 1.012, label = "cluster", angle = 90,
           hjust = 0, fontface = "bold", size = BASE / .pt, family = FONT) +
  scale_fill_gradientn(colours = VIR, limits = c(-2, 2), oob = scales::squish,
                       name = "Row z-score", breaks = c(-2, -1, 0, 1, 2),
                       guide = guide_colourbar(title.position = "top")) +
  scale_x_continuous(limits = c(XB_MIN, XB_MAX), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, nB * 1.105), expand = c(0, 0)) +
  labs(title = paste0("All significant bacterial KOs (n = ", nB, ")"),
       x = NULL, y = NULL) +
  theme_heat() +
  # the panel is now wide enough for a horizontal key underneath, which keeps
  # the grob from growing sideways past the column the layout gives it
  theme(axis.text = element_blank(),
        legend.position = "bottom", legend.justification = "left",
        legend.direction = "horizontal",
        legend.key.width = unit(1.0, "cm"), legend.key.height = unit(0.18, "cm"))
gB <- set_panel_mm(pB, (XB_MAX - XB_MIN) * B_COL_MM, B_H_MM)

# =============================================================================
# C - all fungal KOs with padj < 0.05, tiles + text columns in ONE plot
# =============================================================================
wC <- read.delim(file.path(D, "de_heatmap_fungi_z_wide.tsv"), row.names = 1,
                 check.names = FALSE)[, GROUPS]
hcC <- hclust(dist(as.matrix(wC)), method = "complete")
ordC <- rownames(wC)[hcC$order]; nC <- length(ordC)

DESC <- c(K03696 = "clpC, ATP-dependent Clp protease subunit",
          K00129 = "ALDH3, aldehyde dehydrogenase (NAD(P)+)",
          K15255 = "PIF1, ATP-dependent DNA helicase",
          K00128 = "ALDH, aldehyde dehydrogenase (NAD+)",
          K03884 = "ND6, NADH-ubiquinone oxidoreductase")
PW   <- c(K00129 = "Glycolysis / Gluconeogenesis",
          K00128 = "Glycolysis / Gluconeogenesis",
          K15255 = "DNA replication / repair",
          K03884 = "Oxidative phosphorylation",
          K03696 = "None (Clp chaperone / protease)")
annC <- data.frame(KO = ordC, desc = DESC[ordC], pathway = PW[ordC],
                   yi = match(ordC, rev(ordC)))
zC <- as.data.frame(as.table(as.matrix(wC))); names(zC) <- c("KO", "group", "z")
zC$xi <- match(as.character(zC$group), GROUPS)
zC$yi <- match(as.character(zC$KO), rev(ordC))

X_DESC <- 7.0; X_PATH <- 15.6; XC_MAX <- 21.8
pC <- ggplot() +
  geom_tile(data = zC, aes(xi, yi, fill = z)) +
  geom_text(data = annC, aes(X_DESC, yi, label = desc), hjust = 0,
            size = BASE / .pt, family = FONT, colour = "black") +
  geom_text(data = annC, aes(X_PATH, yi, label = pathway), hjust = 0,
            size = BASE / .pt, family = FONT, colour = "black") +
  annotate("text", x = 1:6, y = nC + 0.75, label = GROUPS, angle = 90, hjust = 0,
           fontface = "bold", size = BASE / .pt, family = FONT) +
  annotate("text", x = c(X_DESC, X_PATH), y = nC + 0.75,
           label = c("KO description", "KEGG pathway"), hjust = 0,
           fontface = "bold", size = BASE / .pt, family = FONT) +
  scale_fill_gradientn(colours = VIR, limits = c(-2, 2), oob = scales::squish,
                       name = "Row z-score", breaks = c(-2, -1, 0, 1, 2)) +
  scale_x_continuous(limits = c(0.4, XC_MAX), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_len(nC), labels = rev(ordC),
                     limits = c(0.4, nC + 1.6), expand = c(0, 0)) +
  labs(title = "All significant fungal KOs (n = 5)", x = NULL, y = NULL) +
  theme_heat() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = BASE, colour = "black"),
        legend.position = "bottom", legend.justification = "left",
        legend.key.width = unit(1.0, "cm"), legend.key.height = unit(0.18, "cm"))
gC <- set_panel_mm(pC, (XC_MAX - 0.4) * C_COL_MM, (nC + 1.2) * C_ROW_MM)

# =============================================================================
# assemble - measure the grobs, hand patchwork those exact millimetres
# =============================================================================
gA <- add_tag(gA, "A"); gB <- add_tag(gB, "B"); gC <- add_tag(gC, "C")
PAD <- 2

h_A <- gh(gA); h_B <- gh(gB); h_C <- gh(gC)
w_L <- max(gw_full(gA), gw_full(gC)); w_R <- gw_full(gB)
# gw vs gw_full: where they differ, a guide box is overflowing its cell and the
# layout must be sized from the larger number or the legend gets cropped
cat(sprintf("  width mm (gw / gw_full)  A %.0f/%.0f  B %.0f/%.0f  C %.0f/%.0f\n",
            gw(gA), gw_full(gA), gw(gB), gw_full(gB), gw(gC), gw_full(gC)))
h_left <- c(h_A, h_C) + PAD
if (sum(h_left) < h_B + PAD) h_left <- h_left * (h_B + PAD) / sum(h_left)

fig <- wrap_elements(full = pad_right(gA, w_L)) + wrap_elements(full = gB) +
       wrap_elements(full = pad_right(gC, w_L)) +
  plot_layout(design  = "AB\nCB",
              widths  = unit(c(w_L + PAD, w_R + PAD), "mm"),
              heights = unit(h_left, "mm")) +
  plot_annotation(theme = theme(plot.margin = margin(2, 2, 2, 2))) &
  # wrap_elements() inherits the default 5.5 pt margin; left in place it steals
  # ~4 mm from every cell and crops the fixed grob
  theme(plot.margin = margin(0, 0, 0, 0))

W_MM <- w_L + w_R + 2 * PAD + 4
H_MM <- sum(h_left) + 4
save_fig(fig, file.path(OUT, "de_counts_and_heatmaps"), W_MM, H_MM)

cat(sprintf("  A %.0f x %.0f | B %.0f x %.0f | C %.0f x %.0f mm\n",
            gw(gA), h_A, gw(gB), h_B, gw(gC), h_C))
cat(sprintf("  panel B: %d rows over %.0f mm = %.3f mm per KO, %d clusters\n",
            nB, B_H_MM, B_H_MM / nB, ncl))
print(prof)
cat("smallest font on the figure:", BASE, "pt\n")
invisible(dev.off())
if (file.exists(file.path(OUT, "Rplots.pdf"))) unlink(file.path(OUT, "Rplots.pdf"))
