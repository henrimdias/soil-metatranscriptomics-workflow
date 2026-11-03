# =============================================================================
# Taxonomic decomposition: contributor turnover vs within-family change, and
# masking of within-family changes in cellular-machinery transcripts
#
#   A  KEGG maps only: contributor turnover vs within-family change,
#      24 -> 72 h. The six maps the gene units of B-C were derived from are coloured
#      as their unit; no gene-list unit is plotted here (it would duplicate its map)
#   B  within-family change of the four machinery units and of their combined
#      signature (circles), and community-level change of the combined
#      signature (squares); filled = all families, open = Bacillaceae removed
#   C  combined-signature within-family index in each sample (vial), with the
#      pooled (both treatments) mean and 95% CI at each time
#
# The plot carries only panel tags, axis titles and category names.
# Every plotted number is re-derived from the raw table by verify_decomposition_plot.py.
#
# Run from this folder:  python3 prepare_decomposition_plot.py ... && Rscript plot_decomposition_masking.R
# =============================================================================
HERE <- Sys.getenv("FIG_DIR", ".")
source(file.path(HERE, "_style.R"))
D <- Sys.getenv("FIG_DATA", file.path(HERE, "data"))

CORE <- c("Ribosomal proteins", "Chaperones (GroEL/DnaK)", "Sec / SRP", "F-ATP synthase")
# PRGn palette for the four machinery units
UNIT_COL <- c("Ribosomal proteins"      = "#40004B",
              "Chaperones (GroEL/DnaK)" = "#9970AB",
              "Sec / SRP"               = "#1B7837",
              "F-ATP synthase"          = "#A6DBA0",
              "Peptidoglycan"           = "grey62",   # cell-wall units: not part of the signature
              "DAP branch"              = "grey62")
SHORT <- c("Ribosomal proteins" = "Ribosomal proteins", "Chaperones (GroEL/DnaK)" = "Chaperones",
           "Sec / SRP" = "Sec/SRP", "F-ATP synthase" = "F-ATP synthase",
           "Peptidoglycan" = "Peptidoglycan", "DAP branch" = "DAP branch")
tag_theme <- theme(plot.tag = element_text(size = BASE + 3, face = "bold", family = FONT),
                   plot.tag.position = c(0, 1))

# =============================================================================
# A
# =============================================================================
maps  <- read.delim(file.path(D, "decomposition_maps.tsv"))
stopifnot(nrow(maps) == 105, all(grepl("^ko", maps$pathway)))      # KEGG maps only
# the maps the gene units of B-C were derived from, one per unit, coloured as the unit
hl <- read.delim(file.path(D, "decomposition_maps_highlight.tsv"))
stopifnot(all(hl$pathway %in% maps$pathway))
hl$name <- c("Ribosome", "RNA degradation", "Protein export", "Photosynthesis",
             "Peptidoglycan\nbiosynthesis", "Lysine\nbiosynthesis")[match(hl$pathway,
             c("ko03010", "ko03018", "ko03060", "ko00195", "ko00550", "ko00300"))]
# label anchors set by hand (Ribosome and RNA degradation lie 0.02 apart)
POS <- data.frame(pathway = c("ko03010", "ko03018", "ko03060", "ko00195", "ko00550", "ko00300"),
                  lx = c(-0.03, 0.60, 0.60, -0.03, -0.40, -0.25),
                  ly = c(-0.40, -0.47, -0.07, -0.25, -0.52, -0.08),
                  hj = c(1, 1, 1, 1, 0.5, 0.5))
hl <- merge(hl, POS, by = "pathway")
grey_maps <- maps[!maps$pathway %in% hl$pathway, ]
LIM <- 0.62
pA <- ggplot() +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.25) +
  geom_point(data = grey_maps, aes(turnover, within), shape = 16, size = 1.0, colour = "grey82") +
  geom_segment(data = hl, aes(turnover, within, xend = lx, yend = ly),
               linewidth = 0.2, colour = "grey45") +
  geom_point(data = hl, aes(turnover, within, fill = unit), shape = 21, size = 2.1,
             colour = "black", stroke = 0.3) +
  geom_label(data = hl, aes(lx, ly, label = name, hjust = hj), size = BASE / .pt,
             family = FONT, label.size = 0, label.padding = unit(0.5, "pt"), fill = "white",
             lineheight = 0.9) +
  scale_fill_manual(values = UNIT_COL, guide = "none") +
  coord_fixed(xlim = c(-LIM, LIM), ylim = c(-LIM, LIM), expand = FALSE) +
  scale_x_continuous(breaks = seq(-0.6, 0.6, 0.3), labels = function(x) sub("-", "\u2212", x)) +
  scale_y_continuous(breaks = seq(-0.6, 0.6, 0.3), labels = function(x) sub("-", "\u2212", x)) +
  labs(tag = "A", x = "Contributor turnover", y = "Within-family change") +
  theme_fig() + tag_theme

# =============================================================================
# B
# =============================================================================
b <- read.delim(file.path(D, "decomposition_units.tsv"))
COMB <- "Combined signature (4 units)"
YPOS <- c("Ribosomal proteins" = 6, "Chaperones (GroEL/DnaK)" = 5, "Sec / SRP" = 4,
          "F-ATP synthase" = 3, "within" = 1.75, "community" = 0.75)
OFF <- 0.17
w <- b |> filter(component == "within", label %in% c(CORE, COMB)) |>
  mutate(row = ifelse(label == COMB, "within", label),
         y = YPOS[row] + ifelse(variant == "baseline", OFF, -OFF),
         col = ifelse(label == COMB, "black", UNIT_COL[label]))
cm <- b |> filter(component == "total", label == COMB) |>
  mutate(y = YPOS["community"] + ifelse(variant == "baseline", OFF, -OFF))
w$fill  <- ifelse(w$variant == "baseline", w$col, "white")
cm$fill <- ifelse(cm$variant == "baseline", BLUE, "white")
ylabs <- c("Ribosomal proteins", "Chaperones", "Sec/SRP", "F-ATP synthase",
           "Combined signature,\nwithin families", "Combined signature,\ncommunity")
# value column (percent change); P values stay in the legend
XA <- 0.19; XN <- 0.35
tabB <- rbind(
  w |> transmute(x = ifelse(variant == "baseline", XA, XN), y = YPOS[row],
                 txt = sprintf("\u2212%d", round(-100 * rel_change)),
                 col = ifelse(label == COMB, "black", "black"),
                 face = ifelse(label == COMB, "bold", "plain")),
  cm |> transmute(x = ifelse(variant == "baseline", XA, XN), y = YPOS["community"],
                  txt = sprintf("\u2212%d", round(-100 * rel_change)), col = BLUE, face = "bold"))
pB <- ggplot() +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.2, ymax = 2.3, fill = "grey94", colour = NA) +
  geom_text(data = tabB, aes(x, y, label = txt, colour = I(col), fontface = face),
            size = BASE / .pt, family = FONT) +
  annotate("text", x = (XA + XN) / 2, y = 7.05, label = "Change (%)", fontface = "bold",
           size = BASE / .pt, family = FONT) +
  annotate("text", x = c(XA, XN), y = 6.62, label = c("All", "No Bac."),
           size = BASE / .pt, family = FONT) +
  annotate("segment", x = XA - 0.08, xend = XN + 0.08, y = 6.42, yend = 6.42, linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.25) +
  geom_errorbarh(data = w, aes(xmin = lo, xmax = hi, y = y, colour = I(col)),
                 height = 0, linewidth = 0.45) +
  geom_point(data = w, aes(rel_change, y, colour = I(col), fill = I(fill)),
             shape = 21, size = 2, stroke = 0.55) +
  geom_errorbarh(data = cm, aes(xmin = lo, xmax = hi, y = y), colour = BLUE,
                 height = 0, linewidth = 0.45) +
  geom_point(data = cm, aes(rel_change, y, fill = I(fill)), shape = 22, size = 2,
             colour = BLUE, stroke = 0.55) +
  scale_x_continuous(limits = c(-0.6, 0.45), breaks = seq(-0.6, 0, 0.2),
                     labels = function(x) sub("-", "−", x), expand = c(0, 0)) +
  scale_y_continuous(breaks = YPOS, labels = ylabs, limits = c(0.2, 7.3), expand = c(0, 0)) +
  labs(tag = "B", x = "Relative change, 24 to 72 h", y = NULL) +
  theme_fig() + tag_theme +
  theme(axis.text.y = element_text(size = BASE, colour = "black", lineheight = 0.9),
        axis.ticks.y = element_blank(), axis.line.y = element_blank())

# =============================================================================
# C
# =============================================================================
fl <- read.delim(file.path(D, "signature_per_sample.tsv"))
mn <- read.delim(file.path(D, "signature_time_means.tsv"))
stopifnot(nrow(fl) == 18, all(table(fl$time) == 6), all(mn$count == 6))
# treatments side by side within each time; small deterministic spread so no point hides another
fl <- fl |> group_by(time, treatment) |> arrange(within_index_rel, .by_group = TRUE) |>
  mutate(x = time + ifelse(treatment == "Control", -4.5, 4.5) + c(-1.2, 0, 1.2)) |> ungroup()
pC <- ggplot() +
  geom_point(data = fl, aes(x, within_index_rel, colour = treatment), size = 1.5) +
  geom_errorbar(data = mn, aes(x = time, ymin = lo, ymax = hi), width = 3, linewidth = 0.4) +
  geom_segment(data = mn, aes(x = time - 3.2, xend = time + 3.2, y = mean, yend = mean),
               linewidth = 0.9) +
  scale_colour_manual(values = c(Control = GREEN, Glucose = PURPLE), name = NULL) +
  # key for the black bars, drawn under the colour key
  annotate("segment", x = 16, xend = 21, y = 0.735, yend = 0.735, linewidth = 0.9) +
  annotate("segment", x = 18.5, xend = 18.5, y = 0.70, yend = 0.77, linewidth = 0.4) +
  annotate("text", x = 23, y = 0.735, label = "Mean \u00b1 95% CI", hjust = 0,
           size = BASE / .pt, family = FONT) +
  scale_x_continuous(breaks = c(24, 48, 72), limits = c(14, 82)) +
  scale_y_continuous(limits = c(0.7, 1.4), breaks = seq(0.7, 1.4, 0.1)) +
  labs(tag = "C", x = "Time (h)", y = "Combined signature index\n(within families)") +
  theme_fig() + tag_theme +
  theme(legend.position = c(1, 1), legend.justification = c(1, 1),
        legend.key.size = unit(0.28, "cm"), legend.background = element_blank(),
        legend.margin = margin(0, 0, 0, 0)) +
  guides(colour = guide_legend(override.aes = list(size = 1.6)))

# =============================================================================
# assemble
# =============================================================================
# shared key: unit colours (A and B) and symbol meaning (B)
k1 <- data.frame(x = c(0, 2.0, 3.25, 4.2, 5.8, 7.5), y = 1,
                 fill = c(UNIT_COL[CORE], "grey62", "grey82"),
                 shp = c(21, 21, 21, 21, 21, 16),
                 lab = c("Ribosomal proteins", "Chaperones", "Sec/SRP", "F-ATP synthase",
                         "Cell wall", "Other KEGG maps"))
k2 <- data.frame(x = c(0, 2.0, 4.7, 6.1), y = 0,
                 fill = c("black", "white", BLUE, "white"), col = c("black", "black", BLUE, BLUE),
                 shp = c(21, 21, 22, 22),
                 lab = c("Within families, all", "Within families, no Bacillaceae",
                         "Community, all", "Community, no Bacillaceae"))
pKey <- ggplot() +
  geom_point(data = k1[k1$shp == 21, ], aes(x, y, fill = I(fill)), shape = 21, size = 2,
             colour = "black", stroke = 0.3) +
  geom_point(data = k1[k1$shp == 16, ], aes(x, y), shape = 16, size = 1.2, colour = "grey82") +
  geom_text(data = k1, aes(x + 0.12, y, label = lab), hjust = 0, size = BASE / .pt, family = FONT) +
  geom_point(data = k2, aes(x, y, fill = I(fill), colour = I(col), shape = I(shp)), size = 2, stroke = 0.55) +
  geom_text(data = k2, aes(x + 0.12, y, label = lab), hjust = 0, size = BASE / .pt, family = FONT) +
  scale_x_continuous(limits = c(-0.2, 10.6), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.6, 1.6)) +
  theme_void()

top <- pA + pB + pC + plot_layout(widths = c(1, 1.05, 0.72))
fig <- top / pKey + plot_layout(heights = c(1, 0.11)) &
  theme(plot.margin = margin(4, 6, 2, 4))
save_fig(fig, file.path(HERE, "decomposition_masking"), 183, 80)
invisible(dev.off())
