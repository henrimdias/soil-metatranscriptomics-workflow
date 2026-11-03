#!/usr/bin/env python3
"""Data for plot_decomposition_masking.R (taxonomic redistribution masks within-family
changes in cellular machinery).
Sources: A2_factorial_contrasts_all_variants.tsv (KEGG maps; decomposition_contrasts.py),
S1_units_factorial.tsv (units) and S3_flask_index.tsv (per sample; machinery_signature.py).
All contrasts: first vs last time point, both treatments pooled, factorial model on all samples.
Run: python3 prepare_decomposition_plot.py <contrasts_dir> <signature_dir> <out_dir>"""
import sys, pathlib
import numpy as np, pandas as pd
from scipy import stats
V, S, OUT = map(pathlib.Path, sys.argv[1:4])
OUT.mkdir(parents=True, exist_ok=True)
C = "time 24->72 (both trt)"
UNITS = {"Ribosomal proteins": "Ribosomal proteins",
         "Chaperones (GroESL / DnaK-DnaJ-GrpE)": "Chaperones (GroEL/DnaK)",
         "Sec / SRP translocon": "Sec / SRP", "F-type ATP synthase": "F-ATP synthase",
         "Cell wall: peptidoglycan": "Peptidoglycan", "Cell wall: DAP branch": "DAP branch",
         "Signature without DAP branch, core 4 units": "Combined signature (4 units)"}
G = pd.read_csv(V / "A2_factorial_contrasts_all_variants.tsv", sep="\t")
m = G[(G.variant == "baseline") & (G.pathway.str.startswith("ko")) & (G.contrast == C) & (G.component != "total")]   # KEGG maps only (drops the Ntr module row)
m = m.pivot_table(index=["pathway", "function"], columns="component", values="est").reset_index()
m.to_csv(OUT / "decomposition_maps.tsv", sep="\t", index=False)

S1 = pd.read_csv(S / "S1_units_factorial.tsv", sep="\t")
S1 = S1[(S1.contrast == C) & S1.unit.isin(UNITS)].copy()
S1["label"] = S1.unit.map(UNITS)
# Panel A shows KEGG maps only. The maps highlighted are the ones the gene units in B-C
# were derived from (one map per unit; Protein export rather than Bacterial secretion
# system for Sec/SRP because it contains all 12 Sec/SRP genes and little else).
HL = pd.DataFrame({"pathway": ["ko03010", "ko03018", "ko03060", "ko00195", "ko00550", "ko00300"],
                   "unit": ["Ribosomal proteins", "Chaperones (GroEL/DnaK)", "Sec / SRP",
                            "F-ATP synthase", "Peptidoglycan", "DAP branch"]})
assert HL.pathway.isin(m.pathway).all()
HL.merge(m, on="pathway").to_csv(OUT / "decomposition_maps_highlight.tsv", sep="\t", index=False)
S1[["variant", "label", "component", "rel_change", "lo", "hi", "p", "q"]].to_csv(OUT / "decomposition_units.tsv", sep="\t", index=False)

S3 = pd.read_csv(S / "S3_flask_index.tsv", sep="\t")
f = S3[(S3.variant == "baseline") & (S3.unit == "Signature without DAP branch, core 4 units")].copy()
f["time"] = f.group.str[1:].astype(int); f["treatment"] = np.where(f.group.str[0] == "C", "Control", "Glucose")
f[["sample", "treatment", "time", "within_index_rel"]].to_csv(OUT / "signature_per_sample.tsv", sep="\t", index=False)
g = f.groupby("time").within_index_rel.agg(["mean", "std", "count"]).reset_index()
t = stats.t.ppf(.975, g["count"] - 1)
g["lo"] = g["mean"] - t * g["std"] / np.sqrt(g["count"]); g["hi"] = g["mean"] + t * g["std"] / np.sqrt(g["count"])
g.to_csv(OUT / "signature_time_means.tsv", sep="\t", index=False)
print(g.round(3).to_string())
print(S1[S1.label.str.startswith("Combined")].round(3).to_string())
