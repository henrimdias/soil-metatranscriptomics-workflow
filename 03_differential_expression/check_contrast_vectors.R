# Prove every contrast vector in deseq2_ko_by_domain.R equals the intended
# difference of group means under ~ treatment * timepoint.
#
# Method: under this design the fitted mean of group g is X_g %*% beta, where
# X_g is that group's row of the model matrix. So the CORRECT contrast vector for
# "group A minus group B" is simply X_A - X_B -- derived, not assumed. The script
# builds each target difference that way and checks it against the hardcoded
# vector element by element. No DESeq2 needed.

meta <- expand.grid(replicate = 1:3,
                    timepoint = factor(c("24h", "48h", "72h"), levels = c("24h", "48h", "72h")),
                    treatment = factor(c("control", "glucose"), levels = c("control", "glucose")))
X <- model.matrix(~ treatment * timepoint, meta)
cat("model matrix columns (order must match DESeq2 resultsNames):\n")
print(colnames(X))

# one representative row per group
row_of <- function(trt, tp) {
  i <- which(meta$treatment == trt & meta$timepoint == tp)[1]
  X[i, ]
}
g <- list(C24 = row_of("control", "24h"), C48 = row_of("control", "48h"),
          C72 = row_of("control", "72h"), G24 = row_of("glucose", "24h"),
          G48 = row_of("glucose", "48h"), G72 = row_of("glucose", "72h"))

# The vectors exactly as written in deseq2_ko_by_domain.R
HARD <- list(
  glucose_G24_vs_C24   = c(0, 1, 0, 0, 0, 0),
  glucose_G48_vs_C48   = c(0, 1, 0, 0, 1, 0),
  glucose_G72_vs_C72   = c(0, 1, 0, 0, 0, 1),
  time_C48_vs_C24      = c(0, 0, 1, 0, 0, 0),
  time_C72_vs_C24      = c(0, 0, 0, 1, 0, 0),
  time_C72_vs_C48      = c(0, 0, -1, 1, 0, 0),
  time_G48_vs_G24      = c(0, 0, 1, 0, 1, 0),
  time_G72_vs_G24      = c(0, 0, 0, 1, 0, 1),
  time_G72_vs_G48      = c(0, 0, -1, 1, -1, 1),
  interaction_48_vs_24 = c(0, 0, 0, 0, 1, 0),
  interaction_72_vs_24 = c(0, 0, 0, 0, 0, 1),
  interaction_72_vs_48 = c(0, 0, 0, 0, -1, 1)
)

# The same differences derived from the model matrix
TRUTH <- list(
  glucose_G24_vs_C24   = g$G24 - g$C24,
  glucose_G48_vs_C48   = g$G48 - g$C48,
  glucose_G72_vs_C72   = g$G72 - g$C72,
  time_C48_vs_C24      = g$C48 - g$C24,
  time_C72_vs_C24      = g$C72 - g$C24,
  time_C72_vs_C48      = g$C72 - g$C48,
  time_G48_vs_G24      = g$G48 - g$G24,
  time_G72_vs_G24      = g$G72 - g$G24,
  time_G72_vs_G48      = g$G72 - g$G48,
  interaction_48_vs_24 = (g$G48 - g$C48) - (g$G24 - g$C24),
  interaction_72_vs_24 = (g$G72 - g$C72) - (g$G24 - g$C24),
  interaction_72_vs_48 = (g$G72 - g$C72) - (g$G48 - g$C48)
)

cat("\n")
bad <- 0
for (nm in names(HARD)) {
  ok <- isTRUE(all.equal(unname(HARD[[nm]]), unname(TRUTH[[nm]])))
  if (!ok) bad <- bad + 1
  cat(sprintf("%s %-22s hardcoded [%s]  derived [%s]\n",
              if (ok) "OK " else "BAD", nm,
              paste(sprintf("%2d", HARD[[nm]]), collapse = " "),
              paste(sprintf("%2d", TRUTH[[nm]]), collapse = " ")))
}

# Independent second check: give beta arbitrary values, compute every group mean,
# and confirm each contrast recovers the intended arithmetic difference of those
# means. This catches an error that a matching vector pair would hide.
set.seed(1)
beta <- c(5.0, 1.7, -0.9, 2.3, 0.4, -1.1)
mu <- sapply(g, function(r) unname(sum(r * beta)))
m <- function(k) unname(mu[[k]])          # drop names, or c(...) below mangles them
target <- c(
  glucose_G24_vs_C24   = m("G24") - m("C24"),
  glucose_G48_vs_C48   = m("G48") - m("C48"),
  glucose_G72_vs_C72   = m("G72") - m("C72"),
  time_C48_vs_C24      = m("C48") - m("C24"),
  time_C72_vs_C24      = m("C72") - m("C24"),
  time_C72_vs_C48      = m("C72") - m("C48"),
  time_G48_vs_G24      = m("G48") - m("G24"),
  time_G72_vs_G24      = m("G72") - m("G24"),
  time_G72_vs_G48      = m("G72") - m("G48"),
  interaction_48_vs_24 = (m("G48") - m("C48")) - (m("G24") - m("C24")),
  interaction_72_vs_24 = (m("G72") - m("C72")) - (m("G24") - m("C24")),
  interaction_72_vs_48 = (m("G72") - m("C72")) - (m("G48") - m("C48")))

cat("\nnumeric check against arbitrary beta:\n")
for (nm in names(HARD)) {
  got <- sum(HARD[[nm]] * beta)
  ok <- isTRUE(all.equal(got, unname(target[nm])))
  if (!ok) bad <- bad + 1
  cat(sprintf("%s %-22s contrast=%8.4f  intended=%8.4f\n",
              if (ok) "OK " else "BAD", nm, got, unname(target[nm])))
}

cat(sprintf("\n%s\n", if (bad == 0) "ALL 12 CONTRASTS VERIFIED (both checks)"
                      else sprintf("%d FAILURES", bad)))
quit(status = if (bad == 0) 0 else 1)
