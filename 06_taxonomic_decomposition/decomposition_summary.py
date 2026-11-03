#!/usr/bin/env python3
"""Summarise the validation: key conclusions under the factorial model (12 df), with a
Welch 6-vs-6 check for the pooled time effect; leave-one-family-out and assignment
variants side by side. Run after decomposition_contrasts.py: python3 decomposition_summary.py <valdir>"""
import sys, pathlib
import numpy as np, pandas as pd
O = pathlib.Path(sys.argv[1])
G = pd.read_csv(O / "A2_factorial_contrasts_all_variants.tsv", sep="\t")
def get(v, fn, c, comp):
    d = G[(G.variant == v) & (G.function == fn) & (G.contrast == c) & (G.component == comp)]
    return d.iloc[0] if len(d) else None
KEYS = [
 ("T1 ribosomal proteins, 24->72: turnover > 0 and within < 0", "Ribosome", "time 24->72 (both trt)",
  lambda t, w: t.est > 0 and w.est < 0 and t.q < .05 and w.q < .05),
 ("T2 peptidoglycan, 24->72: turnover < 0 and within < 0", "Peptidoglycan biosynthesis", "time 24->72 (both trt)",
  lambda t, w: t.est < 0 and w.est < 0 and t.q < .05 and w.q < .05),
 ("T3 glycolysis, 24->72: within > 0 and larger than turnover", "Glycolysis / Gluconeogenesis", "time 24->72 (both trt)",
  lambda t, w: w.est > 0 and w.q < .05 and abs(w.est) > abs(t.est)),
 ("T4 pentose phosphate, 24->72: both > 0", "Pentose phosphate pathway", "time 24->72 (both trt)",
  lambda t, w: t.est > 0 and w.est > 0 and t.q < .05 and w.q < .05),
 ("T5 N metabolism, 24->72: turnover only", "Nitrogen metabolism", "time 24->72 (both trt)",
  lambda t, w: t.q < .05 and w.q >= .05),
 ("G1 glucose 24 h acts within families: glycolysis", "Glycolysis / Gluconeogenesis", "glucose vs control 24 h",
  lambda t, w: w.est > 0 and w.q < .05),
 ("G2 glucose 24 h within: N metabolism", "Nitrogen metabolism", "glucose vs control 24 h",
  lambda t, w: w.est > 0 and w.q < .05),
 ("G3 transience within families: glycolysis (course 24->72)", "Glycolysis / Gluconeogenesis",
  "course (G24->G72)-(C24->C72)", lambda t, w: w.est < 0 and w.q < .05),
 ("G4 glucose sustains peptidoglycan within families (course 24->72)", "Peptidoglycan biosynthesis",
  "course (G24->G72)-(C24->C72)", lambda t, w: w.est > 0 and w.q < .05),
]
rows = []
for v in G.variant.unique():
    for name, fn, c, test in KEYS:
        t, w = get(v, fn, c, "turnover"), get(v, fn, c, "within")
        if t is None:
            continue
        rows.append(dict(variant=v, conclusion=name, turnover=t.est, turnover_q=t.q, within=w.est,
                         within_q=w.q, within_lo=w.lo, within_hi=w.hi, holds=test(t, w)))
    s = G[(G.variant == v) & (G.set_class == "screen")]
    for c in ("time 24->72 (both trt)", "time 24->48 (both trt)"):
        t = s[(s.contrast == c) & (s.component == "turnover")].set_index("function")
        w = s[(s.contrast == c) & (s.component == "within")].set_index("function")
        if len(t) == 0:
            continue
        opp = ((t.q < .05) & (w.q < .05) & (np.sign(t.est) != np.sign(w.est))).sum()
        rows.append(dict(variant=v, conclusion=f"T6 screen maps with opposite significant components, {c}",
                         n=int(opp), n_tested=len(t), holds=opp > 0))
    for c in ("glucose vs control 24 h", "glucose vs control 48 h", "course (G24->G72)-(C24->C72)"):
        t = s[(s.contrast == c) & (s.component == "turnover")]; w = s[(s.contrast == c) & (s.component == "within")]
        if len(t) == 0:
            continue
        rows.append(dict(variant=v, conclusion=f"G5 screen maps significant, {c} (within / turnover)",
                         n=int((w.q < .05).sum()), n_tested=int((t.q < .05).sum()), holds=(w.q < .05).sum() > 0))
K = pd.DataFrame(rows)
K.to_csv(O / "key_conclusions_factorial.tsv", sep="\t", index=False)
piv = K.pivot_table(index="conclusion", columns="variant", values="holds", aggfunc="first")
piv = piv.map(lambda x: "yes" if x is True else ("no" if x is False else ""))
cols = ["baseline", "drop_Bacillaceae"] + [c for c in piv.columns if c.startswith("drop_") and c != "drop_Bacillaceae"] + \
       ["no_placeholder", "order_level", "misassign_10pct", "misassign_20pct", "misassign_30pct"]
piv[[c for c in cols if c in piv.columns]].to_csv(O / "key_conclusions_matrix.tsv", sep="\t")
print(piv[[c for c in cols if c in piv.columns]].to_string())
print(K[K.conclusion.str.startswith(("T6", "G5"))].pivot_table(index="variant", columns="conclusion", values=["n", "n_tested"], aggfunc="first").to_string())
print(K[K.variant.isin(["baseline", "drop_Bacillaceae"]) & K.conclusion.str.match(r"^[TG]\d ")]
      [["variant", "conclusion", "turnover", "turnover_q", "within", "within_q"]].round(3).to_string())
