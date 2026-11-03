#!/usr/bin/env python3
"""Taxonomic decomposition of functional transcript changes: factorial contrasts and robustness checks.

A  Robustness of the decomposition
   A1 exactness and reference point: base-period (Laspeyres), end-period (Paasche) and
      midpoint (Shapley) splits, with the interaction term dw*dphi reported on its own
   A2 inference on the experimental units: each midpoint component equals the difference
      between group means of a per-sample index, so it is tested with a Welch t-test on
      the 3 + 3 independent flasks (no bootstrap). Course differences
      (G72-G24)-(C72-C24) use common weights over the four groups and an OLS contrast
      on the 12 flasks (8 df).
B  Leave-one-family-out: Bacillaceae removed, then each of the 10 largest families
C  Taxonomic assignment: placeholder ('unclassified ...') families removed; order-level
   aggregation; simulated misassignment (a fraction e of every family's reads moved to
   other families in proportion to their transcriptome share, e = 0.1-0.3)
D  Every test treats the 18 flasks as independent units (destructive sampling); nothing
   pairs replicate numbers across time points.
5  Bacillaceae share of bacterial mRNA vs its relative abundance in 16S amplicons.

Run: python3 decomposition_contrasts.py <long.parquet> <KO1-4.txt> <modules_definition.tsv> <outdir>
"""
import os  # paths are read from environment variables (see README)
import sys, pathlib
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

LONG, KO14, MODS, OUT = sys.argv[1:5]
OUT = pathlib.Path(OUT); OUT.mkdir(parents=True, exist_ok=True)
U = pathlib.Path(os.path.join(os.environ.get('MTX_ROOT', '.'), ''))
GROUPS = ["C24", "G24", "C48", "G48", "C72", "G72"]
SAMPLES = [f"{t}_{h}h_{r}" for t in "CG" for h in (24, 48, 72) for r in (1, 2, 3)]
grp = np.array([f"{s[0]}{s.split('_')[1][:-1]}" for s in SAMPLES])
gi = {g: np.where(grp == g)[0] for g in GROUPS}
PAIRS = [("C24", "C48"), ("C24", "C72"), ("G24", "G48"), ("G24", "G72"),
         ("C24", "G24"), ("C48", "G48"), ("C72", "G72")]
COURSES = [("48", ("G24", "G48", "C24", "C48")), ("72", ("G24", "G72", "C24", "C72"))]

# ------------------------------------------------------------------ sets
h = pd.read_csv(KO14, sep="\t")
h["pid"] = h.Pathway.str.extract(r"\[PATH:(ko\d+)\]")
maps = h.dropna(subset=["pid"]).groupby("pid").KO.apply(set)
mname = (h.dropna(subset=["pid"]).drop_duplicates("pid").set_index("pid").Pathway
           .str.replace(r" \[PATH:.*", "", regex=True))
mods = pd.read_csv(MODS, sep="\t")
PRE = ["ko00010", "ko00030", "ko00020", "ko00620", "ko03010", "ko00290", "ko00220",
       "ko00550", "ko00910"]

L = pd.read_parquet(LONG)
fam_dom = {}
for d in ("bacteria", "fungi"):
    for f in pd.read_csv(U / f"inputs/counts_family_{d}.tsv", sep="\t").family:
        fam_dom[f] = d.capitalize()
L["domain"] = L.family.map(fam_dom)
univ = pd.read_csv(U / "results/ko_matrices/counts_ko_bacteria.tsv",
                   sep="\t", index_col=0).index
B = L[(L.domain == "Bacteria") & L.KO.isin(univ)]
CUBE = B.pivot_table(index=["family", "KO"], columns="sample_id", values="count",
                     aggfunc="sum", fill_value=0)[SAMPLES]
FAMS = CUBE.index.get_level_values(0).unique()
WTOT = CUBE.groupby(level="family").sum().reindex(FAMS)            # family x sample

SETS = {}
for pid, kos in maps.items():
    if int(pid[2:]) < 4000 and len(kos & set(univ)) >= 10:
        SETS[(pid, mname[pid], "pre" if pid in PRE else "screen")] = kos
SETS[("module", "Ntr regulon (N limitation)", "pre")] = set(mods[mods.module.str.startswith("Ntr", na=False)].KO)
kolev = CUBE.index.get_level_values(1)
FSET = {}
for key, kos in SETS.items():
    F = CUBE[kolev.isin(kos)].groupby(level="family").sum().reindex(FAMS, fill_value=0)
    if F.sum().min() >= 200:
        FSET[key] = F
print(f"{len(FAMS)} bacterial families, {len(FSET)} function sets "
      f"({sum(k[2] == 'pre' for k in FSET)} pre-declared)", flush=True)


# ------------------------------------------------------------------ core
def shares(Wm, Fm):
    """Wm, Fm: family x sample arrays -> w (family share of transcriptome), phi."""
    w = Wm / Wm.sum(0, keepdims=True)
    phi = np.divide(Fm, Wm, out=np.zeros_like(Fm, dtype=float), where=Wm > 0)
    return w, phi


def welch(x, y):
    """difference mean(y)-mean(x), Welch 95% CI and p"""
    d = y.mean() - x.mean()
    vx, vy = x.var(ddof=1) / len(x), y.var(ddof=1) / len(y)
    se = np.sqrt(vx + vy)
    if se == 0:
        return d, d, d, 1.0
    df = (vx + vy) ** 2 / (vx ** 2 / (len(x) - 1) + vy ** 2 / (len(y) - 1))
    t = stats.t.ppf(0.975, df)
    return d, d - t * se, d + t * se, 2 * stats.t.sf(abs(d / se), df)


def decompose(Wm, Fm):
    w, phi = shares(Wm, Fm)
    out = []
    for a, b in PAIRS:
        wa, wb = w[:, gi[a]].mean(1), w[:, gi[b]].mean(1)
        fa, fb = phi[:, gi[a]].mean(1), phi[:, gi[b]].mean(1)
        Sa, Sb = (wa * fa).sum(), (wb * fb).sum()
        dw, dphi = wb - wa, fb - fa
        T_L, W_L, I = (dw * fa).sum(), (wa * dphi).sum(), (dw * dphi).sum()
        T_P, W_P = (dw * fb).sum(), (wb * dphi).sum()
        T_M, W_M = T_L + I / 2, W_L + I / 2
        # per-flask indices for the midpoint components
        wbar, fbar = (wa + wb) / 2, (fa + fb) / 2
        Iw = (wbar[:, None] * phi).sum(0)        # within index
        Jt = (w * fbar[:, None]).sum(0)          # turnover index
        tw = welch(Iw[gi[a]], Iw[gi[b]]); tt = welch(Jt[gi[a]], Jt[gi[b]])
        out.append(dict(from_group=a, to_group=b, S_from=Sa, S_to=Sb,
                        rel_total=(Sb - Sa) / Sa, exact_err=abs((Sb - Sa) - (T_M + W_M)),
                        turnover_L=T_L / Sa, within_L=W_L / Sa, interaction=I / Sa,
                        turnover_P=T_P / Sa, within_P=W_P / Sa,
                        turnover=T_M / Sa, within=W_M / Sa,
                        turnover_lo=tt[1] / Sa, turnover_hi=tt[2] / Sa, turnover_p=tt[3],
                        within_lo=tw[1] / Sa, within_hi=tw[2] / Sa, within_p=tw[3]))
    for tag, (g1, g2, c1, c2) in COURSES:
        idx = np.concatenate([gi[g] for g in (g1, g2, c1, c2)])
        wbar = np.mean([w[:, gi[g]].mean(1) for g in (g1, g2, c1, c2)], 0)
        fbar = np.mean([phi[:, gi[g]].mean(1) for g in (g1, g2, c1, c2)], 0)
        Sbar = (wbar * fbar).sum()
        for comp, idxv in (("within", (wbar[:, None] * phi).sum(0)),
                           ("turnover", (w * fbar[:, None]).sum(0))):
            m = {g: idxv[gi[g]] for g in (g1, g2, c1, c2)}
            est = (m[g2].mean() - m[g1].mean()) - (m[c2].mean() - m[c1].mean())
            s2 = sum(((v - v.mean()) ** 2).sum() for v in m.values()) / 8
            se = np.sqrt(s2 * 4 / 3)
            t = stats.t.ppf(0.975, 8)
            p = 2 * stats.t.sf(abs(est / se), 8) if se > 0 else 1.0
            out.append(dict(from_group=f"course{tag}", to_group=f"(G24->G{tag})-(C24->C{tag})",
                            S_from=Sbar, **{comp: est / Sbar, f"{comp}_lo": (est - t * se) / Sbar,
                                            f"{comp}_hi": (est + t * se) / Sbar, f"{comp}_p": p}))
    df = pd.DataFrame(out)
    df = df.groupby(["from_group", "to_group"], sort=False).first().reset_index()
    return df


def run_all(Wdf, Fdict, label, only_pre=False):
    rows = []
    fams = Wdf.index
    for (pid, name, cls), F in Fdict.items():
        if only_pre and cls != "pre":
            continue
        d = decompose(Wdf.to_numpy(float), F.reindex(fams, fill_value=0).to_numpy(float))
        d.insert(0, "set_class", cls); d.insert(0, "function", name); d.insert(0, "pathway", pid)
        rows.append(d)
    R = pd.concat(rows, ignore_index=True)
    R.insert(0, "variant", label)
    for comp in ("turnover", "within"):
        R[f"{comp}_q"] = np.nan
        for (cls, a, b), idx in R.groupby([R.set_class, R.from_group, R.to_group]).groups.items():
            R.loc[idx, f"{comp}_q"] = multipletests(R.loc[idx, f"{comp}_p"].fillna(1), method="fdr_bh")[1]
    return R


# ------------------------------------------------------------------ variants
results = []
results.append(run_all(WTOT, FSET, "baseline"))
print("baseline done", flush=True)

# B  leave-one-family-out
share_all = (WTOT / WTOT.sum()).mean(1).sort_values(ascending=False)
top10 = list(share_all.index[:10])
pd.DataFrame({"family": top10, "mean_share": share_all[top10].values}).to_csv(
    OUT / "B_top10_families.tsv", sep="\t", index=False)
for fam in top10:
    keep = [f for f in FAMS if f != fam]
    lab = f"drop_{fam}"
    results.append(run_all(WTOT.loc[keep], {k: v.loc[keep] for k, v in FSET.items()}, lab,
                           only_pre=(fam != "Bacillaceae")))
    print(lab, "done", flush=True)

# C1 placeholder families removed
keep = [f for f in FAMS if not f.lower().startswith("unclassified")]
ph_share = 1 - WTOT.loc[keep].sum() / WTOT.sum()
print(f"placeholder families: {len(FAMS) - len(keep)}, {ph_share.min():.1%}-{ph_share.max():.1%} of assignments")
results.append(run_all(WTOT.loc[keep], {k: v.loc[keep] for k, v in FSET.items()}, "no_placeholder"))

# C2 order-level aggregation (family -> order from the Emu 16S table; placeholders carry the order in the name)
emu = pd.read_csv(U / "inputs/amplicon_16S/relative_emu-combined-genus_16S_BACTERIA.tsv", sep="\t")
f2o = emu.dropna(subset=["family", "order"]).drop_duplicates("family").set_index("family").order.to_dict()
def to_order(f):
    if f in f2o:
        return f2o[f]
    if f.startswith("unclassified ") and f.endswith(" family"):
        return f[len("unclassified "):-len(" family")]
    return f"[{f}]"                       # unmapped family stays its own unit
omap = pd.Series({f: to_order(f) for f in FAMS})
mapped = (~omap.str.startswith("[")).mean()
mapped_w = WTOT.loc[omap.index[~omap.str.startswith("[")]].sum().sum() / WTOT.sum().sum()
print(f"order mapping: {mapped:.1%} of families, {mapped_w:.1%} of assignments")
Wo = WTOT.groupby(omap).sum()
results.append(run_all(Wo, {k: v.groupby(omap).sum() for k, v in FSET.items()}, "order_level"))

# C3 simulated misassignment
for e in (0.1, 0.2, 0.3):
    w = WTOT / WTOT.sum()
    Fe = {k: (1 - e) * v + e * v.sum() * w for k, v in FSET.items()}
    results.append(run_all(WTOT, Fe, f"misassign_{int(e * 100)}pct"))
    print(f"misassign {e} done", flush=True)

R = pd.concat(results, ignore_index=True)
R.to_csv(OUT / "validation_all_variants.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ A1 reference-point check (baseline)
b = R[(R.variant == "baseline") & ~R.from_group.str.startswith("course")].copy()
ref = []
for comp in ("turnover", "within"):
    for alt in ("L", "P"):
        same_sign = (np.sign(b[comp]) == np.sign(b[f"{comp}_{alt}"])) | (b[comp].abs() < 0.01)
        ref.append(dict(component=comp, reference={"L": "base period (Laspeyres)", "P": "end period (Paasche)"}[alt],
                        n=len(b), same_sign_pct=100 * same_sign.mean(),
                        median_abs_diff=(b[comp] - b[f"{comp}_{alt}"]).abs().median()))
dom = lambda T, W: np.where(T.abs() > W.abs(), "turnover", "within")
for alt in ("L", "P"):
    agree = (dom(b.turnover, b.within) == dom(b[f"turnover_{alt}"], b[f"within_{alt}"]))
    ref.append(dict(component="dominant component", reference=alt, n=len(b), same_sign_pct=100 * agree.mean(),
                    median_abs_diff=np.nan))
ref.append(dict(component="|interaction| / (|turnover|+|within|)", reference="midpoint", n=len(b),
                same_sign_pct=np.nan,
                median_abs_diff=(b.interaction.abs() / (b.turnover.abs() + b.within.abs())).median()))
pd.DataFrame(ref).to_csv(OUT / "A1_reference_point.tsv", sep="\t", index=False)
print("exactness: max |dS - (T+W)| =", b.exact_err.max())

# ------------------------------------------------------------------ key conclusions per variant
def row(R, v, fn, a, b_):
    d = R[(R.variant == v) & (R.function == fn) & (R.from_group == a) & (R.to_group == b_)]
    return d.iloc[0] if len(d) else None
KEYS = [
    ("K1 glucose 24 h acts within families: glycolysis", "Glycolysis / Gluconeogenesis", "C24", "G24",
     lambda r: r.within > 0 and r.within_q < 0.05 and abs(r.turnover) < abs(r.within)),
    ("K1b glucose 24 h within: pyruvate", "Pyruvate metabolism", "C24", "G24",
     lambda r: r.within > 0 and r.within_q < 0.05),
    ("K2 transience within families: glycolysis (G-C course 24->72)", "Glycolysis / Gluconeogenesis", "course72",
     "(G24->G72)-(C24->C72)", lambda r: r.within < 0 and r.within_p < 0.05),
    ("K3 ribosomal proteins, control 24->72: turnover > 0, within < 0", "Ribosome", "C24", "C72",
     lambda r: r.turnover > 0 and r.within < 0 and r.turnover_q < 0.05 and r.within_q < 0.05),
    ("K3b ribosomal proteins, glucose 24->72: turnover > 0, within < 0", "Ribosome", "G24", "G72",
     lambda r: r.turnover > 0 and r.within < 0 and r.turnover_q < 0.05 and r.within_q < 0.05),
    ("K4 peptidoglycan, control 24->72: both components negative", "Peptidoglycan biosynthesis", "C24", "C72",
     lambda r: r.turnover < 0 and r.within < 0 and r.turnover_q < 0.05 and r.within_q < 0.05),
    ("K5 peptidoglycan: glucose sustains within (G-C course 24->72 > 0)", "Peptidoglycan biosynthesis", "course72",
     "(G24->G72)-(C24->C72)", lambda r: r.within > 0 and r.within_p < 0.05),
]
kc = []
for v in R.variant.unique():
    for name, fn, a, b_, test in KEYS:
        r = row(R, v, fn, a, b_)
        if r is None:
            continue
        kc.append(dict(variant=v, conclusion=name, turnover=r.turnover, within=r.within,
                       turnover_q=r.get("turnover_q"), within_q=r.get("within_q"),
                       within_p=r.within_p, holds=bool(test(r))))
    s = R[(R.variant == v) & (R.set_class == "screen")]
    for a, b_ in (("C24", "C72"), ("G24", "G72")):
        d = s[(s.from_group == a) & (s.to_group == b_)]
        opp = ((d.turnover_q < 0.05) & (d.within_q < 0.05) & (np.sign(d.turnover) != np.sign(d.within))).sum()
        kc.append(dict(variant=v, conclusion=f"K6 maps with opposite significant components {a}->{b_}",
                       turnover=np.nan, within=np.nan, holds=bool(opp > 0), n_maps=int(opp), n_tested=len(d)))
    d = s[(s.from_group == "C24") & (s.to_group == "G24")]
    kc.append(dict(variant=v, conclusion="K7 glucose 24 h: maps with significant within vs turnover",
                   n_maps=int((d.within_q < 0.05).sum()), n_tested=int((d.turnover_q < 0.05).sum()),
                   holds=bool((d.within_q < 0.05).sum() > (d.turnover_q < 0.05).sum())))
K = pd.DataFrame(kc)
K.to_csv(OUT / "key_conclusions_by_variant.tsv", sep="\t", index=False)


# ------------------------------------------------------------------ A2b factorial model on all 18 flasks
# Common weights (mean over the six groups), so every component is a linear contrast
# of six independent cell means; pooled residual variance, 12 df. Equal-variance
# assumption; the Welch results above are the conservative bound.
CON = {"time 24->72 (both trt)": {"C72": .5, "G72": .5, "C24": -.5, "G24": -.5},
       "time 24->48 (both trt)": {"C48": .5, "G48": .5, "C24": -.5, "G24": -.5},
       "control 24->72": {"C72": 1, "C24": -1}, "glucose 24->72": {"G72": 1, "G24": -1},
       "glucose vs control 24 h": {"G24": 1, "C24": -1},
       "glucose vs control 48 h": {"G48": 1, "C48": -1},
       "glucose vs control 72 h": {"G72": 1, "C72": -1},
       "course (G24->G72)-(C24->C72)": {"G72": 1, "G24": -1, "C72": -1, "C24": 1},
       "course (G24->G48)-(C24->C48)": {"G48": 1, "G24": -1, "C48": -1, "C24": 1}}


def glm(Wm, Fm):
    w, phi = shares(Wm, Fm)
    wbar = np.mean([w[:, gi[g]].mean(1) for g in GROUPS], 0)
    fbar = np.mean([phi[:, gi[g]].mean(1) for g in GROUPS], 0)
    Sbar = (wbar * fbar).sum()
    S = (w * phi).sum(0)
    out = []
    for comp, v in (("within", (wbar[:, None] * phi).sum(0)), ("turnover", (w * fbar[:, None]).sum(0)),
                    ("total", S)):
        cm = {g: v[gi[g]].mean() for g in GROUPS}
        s2 = sum(((v[gi[g]] - cm[g]) ** 2).sum() for g in GROUPS) / 12
        for cname, c in CON.items():
            est = sum(k * cm[g] for g, k in c.items())
            se = np.sqrt(s2 * sum(k * k for k in c.values()) / 3)
            t = stats.t.ppf(0.975, 12)
            out.append(dict(contrast=cname, component=comp, est=est / Sbar, lo=(est - t * se) / Sbar,
                            hi=(est + t * se) / Sbar, p=2 * stats.t.sf(abs(est / se), 12) if se > 0 else 1.0))
    return pd.DataFrame(out)


def glm_all(Wdf, Fdict, label, only_pre=False):
    rows = []
    for (pid, name, cls), F in Fdict.items():
        if only_pre and cls != "pre":
            continue
        d = glm(Wdf.to_numpy(float), F.reindex(Wdf.index, fill_value=0).to_numpy(float))
        d.insert(0, "set_class", cls); d.insert(0, "function", name); d.insert(0, "pathway", pid)
        rows.append(d)
    G = pd.concat(rows, ignore_index=True); G.insert(0, "variant", label)
    G["q"] = np.nan
    for _, idx in G.groupby([G.set_class, G.contrast, G.component]).groups.items():
        G.loc[idx, "q"] = multipletests(G.loc[idx, "p"], method="fdr_bh")[1]
    return G


VAR = {"baseline": (WTOT, FSET, False)}
for fam in top10:
    keep_ = [f for f in FAMS if f != fam]
    VAR[f"drop_{fam}"] = (WTOT.loc[keep_], {k: v.loc[keep_] for k, v in FSET.items()}, fam != "Bacillaceae")
keep_ = [f for f in FAMS if not f.lower().startswith("unclassified")]
VAR["no_placeholder"] = (WTOT.loc[keep_], {k: v.loc[keep_] for k, v in FSET.items()}, False)
VAR["order_level"] = (Wo, {k: v.groupby(omap).sum() for k, v in FSET.items()}, False)
for e in (0.1, 0.2, 0.3):
    w_ = WTOT / WTOT.sum()
    VAR[f"misassign_{int(e * 100)}pct"] = (WTOT, {k: (1 - e) * v + e * v.sum() * w_ for k, v in FSET.items()}, False)
G = pd.concat([glm_all(*v[:2], k, v[2]) for k, v in VAR.items()], ignore_index=True)
G.to_csv(OUT / "A2_factorial_contrasts_all_variants.tsv", sep="\t", index=False)
print("factorial model done", flush=True)

# ------------------------------------------------------------------ 5: Bacillaceae mRNA share vs 16S
meta = pd.read_csv(U / "inputs/amplicon_16S/amplicon_sample_metadata.tsv", sep="\t")
meta["group"] = np.where(meta.C_Pulse.str.contains("no"), "C", "G") + meta.Timepoint_h.astype(str)
bac16 = emu[emu.family == "Bacillaceae"][meta.Barcode].sum() / emu[meta.Barcode].sum()
mr = (WTOT.loc["Bacillaceae"] / WTOT.sum())
five = pd.DataFrame({
    "group": GROUPS,
    "mRNA_share_mean": [mr.iloc[gi[g]].mean() for g in GROUPS],
    "mRNA_share_sd": [mr.iloc[gi[g]].std() for g in GROUPS],
    "16S_rel_abund_mean": [bac16[meta.Barcode[meta.group == g]].mean() for g in GROUPS],
    "16S_rel_abund_sd": [bac16[meta.Barcode[meta.group == g]].std() for g in GROUPS],
    "n_16S": [int((meta.group == g).sum()) for g in GROUPS]})
five["mRNA_to_16S_ratio"] = five.mRNA_share_mean / five["16S_rel_abund_mean"]
five.to_csv(OUT / "Q5_Bacillaceae_mRNA_vs_16S.tsv", sep="\t", index=False)
print(five.round(4).to_string(index=False))
print("done")

# ------------------------------------------------------------------ A2c heteroscedasticity check for the time effect
# Welch on 6 vs 6 flasks (24 h vs 72 h, treatments pooled), common weights over the six
# groups; no equal-variance assumption.
wl = []
for vname in ("baseline", "drop_Bacillaceae", "order_level", "misassign_30pct"):
    Wv, Fv, _ = VAR[vname]
    for (pid, name, cls), F in Fv.items():
        if cls != "pre":
            continue
        w, phi = shares(Wv.to_numpy(float), F.reindex(Wv.index, fill_value=0).to_numpy(float))
        wbar = np.mean([w[:, gi[g]].mean(1) for g in GROUPS], 0)
        fbar = np.mean([phi[:, gi[g]].mean(1) for g in GROUPS], 0)
        Sbar = (wbar * fbar).sum()
        a = np.concatenate([gi["C24"], gi["G24"]]); b_ = np.concatenate([gi["C72"], gi["G72"]])
        for comp, v in (("within", (wbar[:, None] * phi).sum(0)), ("turnover", (w * fbar[:, None]).sum(0))):
            d, lo, hi, p = welch(v[a], v[b_])
            wl.append(dict(variant=vname, function=name, component=comp, est=d / Sbar, lo=lo / Sbar,
                           hi=hi / Sbar, p=p))
WL = pd.DataFrame(wl)
WL["q"] = WL.groupby(["variant", "component"]).p.transform(lambda s: multipletests(s, method="fdr_bh")[1])
WL.to_csv(OUT / "A2c_welch_time24v72_pooled.tsv", sep="\t", index=False)
print(WL[WL.variant == "baseline"].round(3).to_string())
