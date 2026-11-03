#!/usr/bin/env python3
"""Cellular-machinery signature (sensitivity and synthesis, NOT validation).

The five functional units were identified by the masking analysis on these same data
(v03). They are re-defined here as explicit gene lists chosen by gene name, not from
KEGG maps, and frozen before this script was run (machinery_signature_definition.tsv,
.sha256, .frozen_at). Testing them on the same 18 flasks is a sensitivity analysis of
the pattern, not an independent confirmation.

Questions
  1  Is the temporal decline consistent across the units, at community level, with and
     without Bacillaceae?                                    -> S1_units_factorial.tsv
  2  Is it consistent across families? (each family's own share of its transcriptome
     devoted to each unit, 24 h vs 72 h)                     -> S2_family_consistency.tsv
  3  Is it consistent across flasks?                         -> S3_flask_index.tsv
  4  Is it carried by many genes or by a few abundant KOs?   -> S4_gene_contributions.tsv,
                                                                S5_leave_top_KOs_out.tsv
Inference treats the 18 flasks as independent units (factorial cell means, 12 df).

Run: python3 machinery_signature.py <long.parquet> <signature.tsv> <outdir>
"""
import os  # paths are read from environment variables (see README)
import sys, pathlib
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

LONG, SIG, OUT = sys.argv[1:4]
OUT = pathlib.Path(OUT); OUT.mkdir(parents=True, exist_ok=True)
U = pathlib.Path(os.path.join(os.environ.get('MTX_ROOT', '.'), ''))
GROUPS = ["C24", "G24", "C48", "G48", "C72", "G72"]
SAMPLES = [f"{t}_{h}h_{r}" for t in "CG" for h in (24, 48, 72) for r in (1, 2, 3)]
grp = np.array([f"{s[0]}{s.split('_')[1][:-1]}" for s in SAMPLES])
gi = {g: np.where(grp == g)[0] for g in GROUPS}
CON = {"time 24->72 (both trt)": {"C72": .5, "G72": .5, "C24": -.5, "G24": -.5},
       "time 24->48 (both trt)": {"C48": .5, "G48": .5, "C24": -.5, "G24": -.5},
       "control 24->72": {"C72": 1, "C24": -1}, "glucose 24->72": {"G72": 1, "G24": -1}}
T24 = np.concatenate([gi["C24"], gi["G24"]]); T48 = np.concatenate([gi["C48"], gi["G48"]])
T72 = np.concatenate([gi["C72"], gi["G72"]])

sig = pd.read_csv(SIG, sep="\t")
UNITS = {u: list(d.KO) for u, d in sig.groupby("unit", sort=False)}
UNITS["Signature (all units)"] = list(sig.KO)
UNITS["Signature without DAP branch"] = list(sig[sig.unit != "Cell wall: DAP branch"].KO)
# post hoc, after the gene-level check showed both cell-wall units are carried by 1-3 KOs
UNITS["Signature without DAP branch, core 4 units"] = list(sig[~sig.unit.str.startswith("Cell wall")].KO)
SYM = dict(zip(sig.KO, sig.symbol)); UNIT_OF = dict(zip(sig.KO, sig.unit))

L = pd.read_parquet(LONG)
fam = set(pd.read_csv(U / "inputs/counts_family_bacteria.tsv", sep="\t").family)
univ = pd.read_csv(U / "results/ko_matrices/counts_ko_bacteria.tsv", sep="\t", index_col=0).index
L = L[L.family.isin(fam) & L.KO.isin(univ)]
CUBE = L.pivot_table(index=["family", "KO"], columns="sample_id", values="count", aggfunc="sum", fill_value=0)[SAMPLES]
WTOT = CUBE.groupby(level="family").sum()
kolev = CUBE.index.get_level_values(1)
KOCOUNT = {k: CUBE[kolev == k].droplevel(1) for k in sig.KO if (kolev == k).any()}
print(f"{len(WTOT)} families; signature KOs present: {len(KOCOUNT)}/{len(sig)}", flush=True)


def contrast_stats(v):
    cm = {g: v[gi[g]].mean() for g in GROUPS}
    s2 = sum(((v[gi[g]] - cm[g]) ** 2).sum() for g in GROUPS) / 12
    out = {}
    for name, c in CON.items():
        est = sum(k * cm[g] for g, k in c.items())
        se = np.sqrt(s2 * sum(k * k for k in c.values()) / 3)
        out[name] = (est, se, 2 * stats.t.sf(abs(est / se), 12) if se > 0 else 1.0)
    return out


def components(Wdf, Fdf):
    """per-flask indices with weights common to the six groups"""
    Wm, Fm = Wdf.to_numpy(float), Fdf.reindex(Wdf.index, fill_value=0).to_numpy(float)
    w = Wm / Wm.sum(0, keepdims=True)
    phi = np.divide(Fm, Wm, out=np.zeros_like(Fm), where=Wm > 0)
    wbar = np.mean([w[:, gi[g]].mean(1) for g in GROUPS], 0)
    fbar = np.mean([phi[:, gi[g]].mean(1) for g in GROUPS], 0)
    Sbar = (wbar * fbar).sum()
    return dict(total=(w * phi).sum(0), within=(wbar[:, None] * phi).sum(0),
                turnover=(w * fbar[:, None]).sum(0)), Sbar, wbar, phi


VARIANTS = {"baseline": list(WTOT.index), "drop_Bacillaceae": [f for f in WTOT.index if f != "Bacillaceae"]}

# ---------------------------------------------------------------- 1 units, factorial
r1, r3 = [], []
for v, keep in VARIANTS.items():
    Wv = WTOT.loc[keep]
    for u, kos in UNITS.items():
        F = sum(KOCOUNT[k].reindex(keep, fill_value=0) for k in kos if k in KOCOUNT)
        comp, Sbar, _, _ = components(Wv, F)
        for cname, idx in comp.items():
            for con, (est, se, p) in contrast_stats(idx).items():
                r1.append(dict(variant=v, unit=u, component=cname, contrast=con, rel_change=est / Sbar,
                               lo=(est - stats.t.ppf(.975, 12) * se) / Sbar,
                               hi=(est + stats.t.ppf(.975, 12) * se) / Sbar, p=p))
        if u.startswith("Signature"):
            for j, s in enumerate(SAMPLES):
                r3.append(dict(variant=v, unit=u, sample=s, group=grp[j],
                               within_index_rel=comp["within"][j] / Sbar,
                               total_share=comp["total"][j]))
S1 = pd.DataFrame(r1)
S1["q"] = np.nan
base_units = ~S1.unit.str.startswith("Signature")
for _, idx in S1[base_units].groupby(["variant", "component", "contrast"]).groups.items():
    S1.loc[idx, "q"] = multipletests(S1.loc[idx, "p"], method="fdr_bh")[1]
S1.loc[~base_units, "q"] = S1.loc[~base_units, "p"]          # composite: single pre-set test
S1.to_csv(OUT / "S1_units_factorial.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- 3 flasks
S3 = pd.DataFrame(r3)
sep = []
for (v, u), d in S3.groupby(["variant", "unit"]):
    x = d.set_index("sample").within_index_rel
    a, b48, b72 = x.iloc[T24].values, x.iloc[T48].values, x.iloc[T72].values
    sep.append(dict(variant=v, unit=u,
                    AUC_24_above_72=float(np.mean([[i > j for j in b72] for i in a])),
                    AUC_24_above_48=float(np.mean([[i > j for j in b48] for i in a])),
                    min_24h=a.min(), max_48h=b48.max(), max_72h=b72.max()))
S3.to_csv(OUT / "S3_flask_index.tsv", sep="\t", index=False)
S3s = pd.DataFrame(sep); S3s.to_csv(OUT / "S3_flask_separation.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- 2 families
share = (WTOT / WTOT.sum()).mean(1)
big = share[share >= 0.005].index
r2 = []
for u, kos in UNITS.items():
    F = sum(KOCOUNT[k].reindex(WTOT.index, fill_value=0) for k in kos if k in KOCOUNT)
    for f in big:
        cnt = F.loc[f].to_numpy(float); tot = WTOT.loc[f].to_numpy(float)
        if cnt[T24].mean() < 20 or cnt[T72].mean() < 20:
            continue
        phi = cnt / tot
        a, b = np.log2(phi[T24] + 1e-9), np.log2(phi[T72] + 1e-9)
        t = stats.ttest_ind(b, a, equal_var=False)
        r2.append(dict(unit=u, family=f, mean_share=share[f], phi_24=phi[T24].mean(), phi_72=phi[T72].mean(),
                       log2FC_72v24=b.mean() - a.mean(), p=t.pvalue,
                       counts_24=cnt[T24].mean(), counts_72=cnt[T72].mean()))
S2 = pd.DataFrame(r2)
S2["q"] = S2.groupby("unit").p.transform(lambda s: multipletests(s, method="fdr_bh")[1])
S2.to_csv(OUT / "S2_family_consistency.tsv", sep="\t", index=False)
S2s = S2.groupby("unit").apply(lambda d: pd.Series({
    "n_families": len(d), "n_down": int((d.log2FC_72v24 < 0).sum()),
    "n_down_q05": int(((d.log2FC_72v24 < 0) & (d.q < .05)).sum()),
    "n_up_q05": int(((d.log2FC_72v24 > 0) & (d.q < .05)).sum()),
    "sign_test_p": stats.binomtest(int((d.log2FC_72v24 < 0).sum()), len(d), 0.5).pvalue,
    "median_log2FC": d.log2FC_72v24.median(),
    "n_down_excl_Bacillaceae": int(((d.log2FC_72v24 < 0) & (d.family != "Bacillaceae")).sum()),
    "n_excl_Bacillaceae": int((d.family != "Bacillaceae").sum())}), include_groups=False)
S2s.to_csv(OUT / "S2_family_consistency_summary.tsv", sep="\t")

# ---------------------------------------------------------------- 4 genes
r4, r5 = [], []
for v, keep in VARIANTS.items():
    Wv = WTOT.loc[keep]
    for u, kos in UNITS.items():
        if u.startswith("Signature without DAP branch") and "core" not in u:
            continue
        present = [k for k in kos if k in KOCOUNT]
        Fu = sum(KOCOUNT[k].reindex(keep, fill_value=0) for k in present)
        compu, Sbar, wbar, _ = components(Wv, Fu)
        unit_est = contrast_stats(compu["within"])["time 24->72 (both trt)"][0]
        for k in present:
            ck, _, _, _ = components(Wv, KOCOUNT[k].reindex(keep, fill_value=0))
            # same unit-level weights so KO contributions add up to the unit's within change
            Wm = Wv.to_numpy(float)
            phik = np.divide(KOCOUNT[k].reindex(keep, fill_value=0).to_numpy(float), Wm,
                             out=np.zeros(Wm.shape), where=Wm > 0)
            idx = (wbar[:, None] * phik).sum(0)
            est, se, p = contrast_stats(idx)["time 24->72 (both trt)"]
            r4.append(dict(variant=v, unit=u, KO=k, symbol=SYM[k], member_unit=UNIT_OF[k],
                           share_of_unit=float(idx.mean() / compu["within"].mean()),
                           within_contrib=est / Sbar, ko_rel_change=est / idx.mean() if idx.mean() > 0 else np.nan,
                           p=p, frac_of_unit_change=est / unit_est if unit_est != 0 else np.nan))
        # leave the 1 / 3 / 5 largest contributors out
        d = pd.DataFrame([x for x in r4 if x["variant"] == v and x["unit"] == u])
        order = d.sort_values("within_contrib").KO.tolist()            # most negative first
        for n in (0, 1, 3, 5):
            ks = [k for k in present if k not in order[:n]]
            if not ks:
                continue
            F = sum(KOCOUNT[k].reindex(keep, fill_value=0) for k in ks)
            comp, Sb, _, _ = components(Wv, F)
            est, se, p = contrast_stats(comp["within"])["time 24->72 (both trt)"]
            r5.append(dict(variant=v, unit=u, n_removed=n, removed=", ".join(SYM[k] for k in order[:n]),
                           within_rel_change=est / Sb, p=p))
S4 = pd.DataFrame(r4)
S4["q"] = S4.groupby(["variant", "unit"]).p.transform(lambda s: multipletests(s, method="fdr_bh")[1])
S4.to_csv(OUT / "S4_gene_contributions.tsv", sep="\t", index=False)
S5 = pd.DataFrame(r5); S5.to_csv(OUT / "S5_leave_top_KOs_out.tsv", sep="\t", index=False)

conc = []
for (v, u), d in S4.groupby(["variant", "unit"]):
    neg = d[d.within_contrib < 0].sort_values("within_contrib")
    tot_neg = neg.within_contrib.sum()
    cum = (neg.within_contrib.cumsum() / tot_neg).to_numpy()
    conc.append(dict(variant=v, unit=u, n_KO=len(d), n_KO_down=len(neg),
                     n_KO_down_q05=int(((d.within_contrib < 0) & (d.q < .05)).sum()),
                     n_KO_up_q05=int(((d.within_contrib > 0) & (d.q < .05)).sum()),
                     KOs_for_50pct=int(np.searchsorted(cum, .5) + 1) if len(cum) else 0,
                     KOs_for_80pct=int(np.searchsorted(cum, .8) + 1) if len(cum) else 0,
                     top3=", ".join(neg.symbol.head(3))))
pd.DataFrame(conc).to_csv(OUT / "S4_gene_concentration_summary.tsv", sep="\t", index=False)

pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
x = S1[(S1.contrast == "time 24->72 (both trt)")].pivot_table(index=["unit"], columns=["variant", "component"], values="rel_change")
print(x.round(3).to_string())
print(S1[(S1.contrast == "time 24->72 (both trt)") & (S1.component == "within")][["variant", "unit", "rel_change", "lo", "hi", "q"]].round(3).to_string())
print(S2s.round(3).to_string())
print(S3s.round(3).to_string())
print(pd.DataFrame(conc).to_string())
print(S5.round(3).to_string())
