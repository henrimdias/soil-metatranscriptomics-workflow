#!/usr/bin/env python3
"""Summary tables of the taxonomic decomposition of bacterial functional transcript changes.

Recomputes every value from the family x KO x sample counts
(long.parquet, built from counts_ko_taxonomy_full.tsv) and writes one TSV per sheet.
Same estimator as decomposition_contrasts.py (factorial block) and verify_decomposition_plot.py:
common weights = mean of the six group means, pooled residual variance, 12 df.

Run: python3 decomposition_tables.py <long.parquet> <KO1-4.txt> <machinery_signature_definition.tsv>
                            <counts_family_bacteria.tsv> <counts_ko_bacteria.tsv> <emu16S.tsv> <outdir>
"""
import sys, pathlib
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

LONG, KO14, SIG, FAMB, KOB, EMU, OUT = sys.argv[1:8]
OUT = pathlib.Path(OUT); OUT.mkdir(parents=True, exist_ok=True)

GROUPS = ["C24", "G24", "C48", "G48", "C72", "G72"]
SAMPLES = [f"{t}_{h}h_{r}" for t in "CG" for h in (24, 48, 72) for r in (1, 2, 3)]
grp = np.array([f"{s[0]}{s.split('_')[1][:-1]}" for s in SAMPLES])
gi = {g: np.where(grp == g)[0] for g in GROUPS}
CON = {  # name: (label, coefficients on cell means)
    "time 24->72 (both trt)": ("24 to 72 h, both treatments", {"C24": -.5, "G24": -.5, "C72": .5, "G72": .5}),
    "time 24->48 (both trt)": ("24 to 48 h, both treatments", {"C24": -.5, "G24": -.5, "C48": .5, "G48": .5}),
    "control 24->72": ("24 to 72 h, control", {"C24": -1, "C72": 1}),
    "glucose 24->72": ("24 to 72 h, glucose", {"G24": -1, "G72": 1}),
    "glucose vs control 24 h": ("Glucose vs control, 24 h", {"G24": 1, "C24": -1}),
    "glucose vs control 48 h": ("Glucose vs control, 48 h", {"G48": 1, "C48": -1}),
    "glucose vs control 72 h": ("Glucose vs control, 72 h", {"G72": 1, "C72": -1}),
    "interaction 24->48": ("Treatment x time, 24 to 48 h", {"G48": 1, "G24": -1, "C48": -1, "C24": 1}),
    "interaction 24->72": ("Treatment x time, 24 to 72 h", {"G72": 1, "G24": -1, "C72": -1, "C24": 1}),
}
TCRIT = stats.t.ppf(0.975, 12)

# ------------------------------------------------------------------ data
fam_b = set(pd.read_csv(FAMB, sep="\t").family)
univ = pd.read_csv(KOB, sep="\t", index_col=0)
L = pd.read_parquet(LONG)
B = L[L.family.isin(fam_b) & L.KO.isin(univ.index)]
CUBE = B.pivot_table(index=["family", "KO"], columns="sample_id", values="count",
                     aggfunc="sum", fill_value=0)[SAMPLES]
assert (CUBE.groupby(level="KO").sum().reindex(univ.index, fill_value=0)[SAMPLES].values
        == univ[SAMPLES].values).all(), "cube does not reproduce the DESeq2 KO matrix"
FAMS = CUBE.index.get_level_values(0).unique()
WTOT = CUBE.groupby(level="family").sum().reindex(FAMS)
kolev = CUBE.index.get_level_values(1)

h = pd.read_csv(KO14, sep="\t")
h["pid"] = h.Pathway.str.extract(r"\[PATH:(ko\d+)\]")
hp = h.dropna(subset=["pid"])
maps = hp.groupby("pid").KO.apply(set)
minfo = hp.drop_duplicates("pid").set_index("pid")
minfo["name"] = minfo.Pathway.str.replace(r" \[PATH:.*", "", regex=True)
kodesc = h.drop_duplicates("KO").set_index("KO").KoDescription

def fcount(kos, W=WTOT):
    return CUBE[kolev.isin(kos)].groupby(level="family").sum().reindex(W.index, fill_value=0)

MAPS = {}
for pid, kos in maps.items():
    k = kos & set(univ.index)
    if int(pid[2:]) < 4000 and len(k) >= 10:
        F = fcount(k)
        if F.sum().min() >= 200:
            MAPS[pid] = (k, F)
assert len(MAPS) == 105, len(MAPS)

sig = pd.read_csv(SIG, sep="\t")
UNAME = {"Ribosomal proteins": "Ribosomal proteins",
         "Chaperones (GroESL / DnaK-DnaJ-GrpE)": "Chaperones (GroESL, DnaK-DnaJ-GrpE)",
         "Sec / SRP translocon": "Sec/SRP translocon", "F-type ATP synthase": "F-type ATP synthase",
         "Cell wall: peptidoglycan": "Peptidoglycan synthesis", "Cell wall: DAP branch": "Diaminopimelate branch"}
sig["Unit"] = sig.unit.map(UNAME)
UNITS = {u: set(d.KO) for u, d in sig.groupby("Unit", sort=False)}
CORE = ["Ribosomal proteins", "Chaperones (GroESL, DnaK-DnaJ-GrpE)", "Sec/SRP translocon", "F-type ATP synthase"]
COMB = "Combined signature (4 units)"
UNITS[COMB] = set().union(*[UNITS[u] for u in CORE])
UORDER = CORE + ["Peptidoglycan synthesis", "Diaminopimelate branch", COMB]
assert len(UNITS[COMB]) == 105

# ------------------------------------------------------------------ estimator
def decomp(Wm, Fm):
    w = Wm / Wm.sum(0, keepdims=True)
    phi = np.divide(Fm, Wm, out=np.zeros_like(Fm, dtype=float), where=Wm > 0)
    wbar = np.mean([w[:, gi[g]].mean(1) for g in GROUPS], 0)
    fbar = np.mean([phi[:, gi[g]].mean(1) for g in GROUPS], 0)
    Sbar = (wbar * fbar).sum()
    idx = {"Community": (w * phi).sum(0), "Within-family": (wbar[:, None] * phi).sum(0),
           "Turnover": (w * fbar[:, None]).sum(0)}
    rows = []
    for comp, v in idx.items():
        cm = {g: v[gi[g]].mean() for g in GROUPS}
        s2 = sum(((v[gi[g]] - cm[g]) ** 2).sum() for g in GROUPS) / 12
        for cname, (lab, c) in CON.items():
            est = sum(k * cm[g] for g, k in c.items())
            se = np.sqrt(s2 * sum(k * k for k in c.values()) / 3)
            rows.append(dict(contrast=cname, contrast_label=lab, component=comp, est=est / Sbar,
                             lo=(est - TCRIT * se) / Sbar, hi=(est + TCRIT * se) / Sbar,
                             p=2 * stats.t.sf(abs(est / se), 12) if se > 0 else 1.0))
    return pd.DataFrame(rows), Sbar, {k: v / Sbar for k, v in idx.items()}

def bh(df, by):
    df["q"] = df.groupby(by).p.transform(lambda s: multipletests(s, method="fdr_bh")[1])
    return df

# ------------------------------------------------------------------ variants
share_all = (WTOT / WTOT.sum()).mean(1).sort_values(ascending=False)
top10 = list(share_all.index[:10])
emu = pd.read_csv(EMU, sep="\t")
f2o = emu.dropna(subset=["family", "order"]).drop_duplicates("family").set_index("family").order.to_dict()
def to_order(f):
    if f in f2o: return f2o[f]
    if f.startswith("unclassified ") and f.endswith(" family"): return f[len("unclassified "):-len(" family")]
    return f"[{f}]"
omap = pd.Series({f: to_order(f) for f in FAMS})

def variant(name):
    """returns (W, fn(F)->F') for a variant"""
    if name == "All families": return WTOT, lambda F: F
    if name.startswith("Without "):
        keep = [f for f in FAMS if f != name[len("Without "):]]
        return WTOT.loc[keep], lambda F: F.loc[keep]
    if name == "Without unclassified placeholders":
        keep = [f for f in FAMS if not f.lower().startswith("unclassified")]
        return WTOT.loc[keep], lambda F: F.loc[keep]
    if name == "Order level":
        return WTOT.groupby(omap).sum(), lambda F: F.groupby(omap).sum()
    if name.startswith("Misassignment"):
        e = int(name.split()[1].rstrip("%")) / 100
        w_ = WTOT / WTOT.sum()
        return WTOT, lambda F: (1 - e) * F + e * F.sum() * w_
    raise ValueError(name)

VARIANTS = (["All families"] + [f"Without {f}" for f in top10] +
            ["Without unclassified placeholders", "Order level",
             "Misassignment 10%", "Misassignment 20%", "Misassignment 30%"])

def run_sets(sets, vname):
    W, tf = variant(vname)
    out, S = [], {}
    for key, kos in sets.items():
        F = tf(fcount(kos))
        d, Sbar, idx = decomp(W.to_numpy(float), F.reindex(W.index, fill_value=0).to_numpy(float))
        d.insert(0, "set", key); out.append(d); S[key] = (Sbar, idx)
    return pd.concat(out, ignore_index=True), S

mapsets = {pid: v[0] for pid, v in MAPS.items()}
UM, US = {}, {}
for v in VARIANTS:
    m, _ = run_sets(mapsets, v)
    UM[v] = bh(m, ["contrast", "component"])
    u, S = run_sets({u: UNITS[u] for u in UORDER}, v)
    six = u[u.set != COMB].copy(); six = bh(six, ["contrast", "component"])
    u = u.merge(six[["set", "contrast", "component", "q"]], how="left")
    UM[v + "|u"] = u; US[v] = S
    print(v, "done", flush=True)

# ------------------------------------------------------------------ S10a maps (All families, 24->72)
def wide(df):
    d = df[df.contrast == "time 24->72 (both trt)"]
    out = None
    for comp in ("Turnover", "Within-family", "Community"):
        x = d[d.component == comp].set_index("set")[["est", "lo", "hi", "p", "q"]]
        x.columns = [f"{comp}|{c}" for c in x.columns]
        out = x if out is None else out.join(x)
    return out

A = wide(UM["All families"])
_, Sm = run_sets(mapsets, "All families")
A["mean_share"] = [Sm[p][0] for p in A.index]
A["n_KO"] = [len(MAPS[p][0]) for p in A.index]
A["min_assign"] = [int(MAPS[p][1].sum().min()) for p in A.index]
A["residual"] = A["Community|est"] - A["Turnover|est"] - A["Within-family|est"]
A["opposing"] = ((A["Turnover|q"] < .05) & (A["Within-family|q"] < .05) &
                 (np.sign(A["Turnover|est"]) != np.sign(A["Within-family|est"])))
# KOs contributing most to the within-family change of each map (share of the 24->72 within-family change)
def top_kos(kos, n=3):
    W = WTOT.to_numpy(float); w = W / W.sum(0, keepdims=True)
    wbar = np.mean([w[:, gi[g]].mean(1) for g in GROUPS], 0)
    c = {}
    for k in kos:
        F = fcount({k}).to_numpy(float)
        phi = np.divide(F, W, out=np.zeros_like(F), where=W > 0)
        iv = (wbar[:, None] * phi).sum(0)
        c[k] = (iv[np.r_[gi["C72"], gi["G72"]]].mean() - iv[np.r_[gi["C24"], gi["G24"]]].mean())
    s = pd.Series(c); tot = s.sum()
    s = s.reindex(s.abs().sort_values(ascending=False).index)[:n]
    sym = lambda k: str(kodesc.get(k, k)).split(";")[0].split(",")[0]
    return "; ".join(f"{k} {sym(k)} ({100 * v / tot:.0f}%)" for k, v in s.items()) if tot != 0 else ""
A["top_KOs"] = [top_kos(MAPS[p][0]) for p in A.index]
A["L1"] = minfo.loc[A.index, "PathwayL1"].values
A["L2"] = minfo.loc[A.index, "PathwayL2"].values
A["name"] = minfo.loc[A.index, "name"].values
A = A.reset_index().rename(columns={"set": "pid"}).sort_values("pid")
A.to_csv(OUT / "S10a.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ S10b gene lists
sig["KoDescription"] = sig.KO.map(kodesc)
sig["in_signature"] = sig.Unit.isin(CORE)
sig["maps"] = sig.KO.map(lambda k: "; ".join(sorted(p for p, (ks, _) in MAPS.items() if k in ks)))
sig[["Unit", "KO", "symbol", "KoDescription", "in_signature", "maps"]].to_csv(OUT / "S10b.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ S10c units, all contrasts, two variants
C = []
for v in ("All families", "Without Bacillaceae"):
    u = UM[v + "|u"].copy(); u.insert(0, "variant", v)
    u["mean_share"] = u.set.map(lambda s: US[v][s][0])
    C.append(u)
C = pd.concat(C, ignore_index=True)
C.to_csv(OUT / "S10c.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ S10d families (Welch, log2 phi, 6 vs 6)
T24, T72 = np.r_[gi["C24"], gi["G24"]], np.r_[gi["C72"], gi["G72"]]
share = (WTOT / WTOT.sum()).mean(1)
big = share[share >= 0.005].index
r = []
for u in UORDER:
    F = fcount(UNITS[u])
    for f in big:
        cnt, tot = F.loc[f].to_numpy(float), WTOT.loc[f].to_numpy(float)
        if cnt[T24].mean() < 20 or cnt[T72].mean() < 20: continue
        phi = cnt / tot
        a, b = np.log2(phi[T24] + 1e-9), np.log2(phi[T72] + 1e-9)
        t = stats.ttest_ind(b, a, equal_var=False)
        r.append(dict(unit=u, family=f, mean_share=share[f], counts_24=cnt[T24].mean(), counts_72=cnt[T72].mean(),
                      phi_24=phi[T24].mean(), phi_72=phi[T72].mean(), log2FC=b.mean() - a.mean(), p=t.pvalue))
D = pd.DataFrame(r); D = bh(D, ["unit"])
D.to_csv(OUT / "S10d.tsv", sep="\t", index=False)
Ds = D.groupby("unit", sort=False).apply(lambda d: pd.Series(dict(
    n=len(d), n_down=int((d.log2FC < 0).sum()), n_down_q=int(((d.log2FC < 0) & (d.q < .05)).sum()),
    n_up_q=int(((d.log2FC > 0) & (d.q < .05)).sum()),
    sign_p=stats.binomtest(int((d.log2FC < 0).sum()), len(d), .5).pvalue,
    n_noBac=int((d.family != "Bacillaceae").sum()),
    n_down_noBac=int(((d.log2FC < 0) & (d.family != "Bacillaceae")).sum()))), include_groups=False)
Ds.to_csv(OUT / "S10d_summary.tsv", sep="\t")

# ------------------------------------------------------------------ S10e vials
E = []
for v in ("All families", "Without Bacillaceae"):
    for u in UORDER:
        Sbar, idx = US[v][u]
        for i, s in enumerate(SAMPLES):
            E.append(dict(variant=v, unit=u, sample=s, treatment="Control" if s[0] == "C" else "Glucose",
                          time=int(s.split("_")[1][:-1]), within=idx["Within-family"][i],
                          turnover=idx["Turnover"][i], community=idx["Community"][i]))
E = pd.DataFrame(E); E.to_csv(OUT / "S10e.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ S10f robustness of the units (24->72, all variants)
F_ = []
for v in VARIANTS:
    u = UM[v + "|u"]; u = u[u.contrast == "time 24->72 (both trt)"].copy(); u.insert(0, "variant", v); F_.append(u)
pd.concat(F_, ignore_index=True).to_csv(OUT / "S10f.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ S10g robustness of the map-level result
OPP = list(A.loc[A.opposing & (A["Turnover|est"] > 0), "pid"])
G = []
for v in VARIANTS:
    m = UM[v]; m = m[m.contrast == "time 24->72 (both trt)"]
    t = m[m.component == "Turnover"].set_index("set"); w = m[m.component == "Within-family"].set_index("set")
    opp = (t.q < .05) & (w.q < .05) & (np.sign(t.est) != np.sign(w.est))
    pos = opp & (t.est > 0)
    row = dict(variant=v, n_maps=len(t), n_opposing=int(opp.sum()), n_masking=int(pos.sum()))
    for p in OPP: row[p] = bool(pos.get(p, False))
    G.append(row)
G = pd.DataFrame(G); G.to_csv(OUT / "S10g.tsv", sep="\t", index=False)
pd.Series(OPP).to_csv(OUT / "S10g_maps.tsv", sep="\t", index=False, header=False)

# ------------------------------------------------------------------ S10h exact pairwise splits (maps, All families)
# base-period, end-period and midpoint (Shapley) splits of the change between group means
w = (WTOT / WTOT.sum()).to_numpy(float)
H = []
for pid, (kos, Fm) in MAPS.items():
    phi = np.divide(Fm.to_numpy(float), WTOT.to_numpy(float), out=np.zeros(Fm.shape), where=WTOT.to_numpy() > 0)
    for a_, b_ in (("24", "72"), ("24", "48")):
        A_ = np.r_[gi["C" + a_], gi["G" + a_]]; B_ = np.r_[gi["C" + b_], gi["G" + b_]]
        wa, wb, fa, fb = w[:, A_].mean(1), w[:, B_].mean(1), phi[:, A_].mean(1), phi[:, B_].mean(1)
        Sa, Sb = (wa * fa).sum(), (wb * fb).sum()
        dw, df_ = wb - wa, fb - fa
        I = (dw * df_).sum()
        H.append(dict(pid=pid, contrast=f"{a_} to {b_} h", S_from=Sa, S_to=Sb, total=(Sb - Sa) / Sa,
                      turnover_base=(dw * fa).sum() / Sa, within_base=(wa * df_).sum() / Sa,
                      turnover_end=(dw * fb).sum() / Sa, within_end=(wb * df_).sum() / Sa,
                      cross=I / Sa, turnover_mid=((dw * fa).sum() + I / 2) / Sa,
                      within_mid=((wa * df_).sum() + I / 2) / Sa))
H = pd.DataFrame(H)
assert np.allclose(H.turnover_mid + H.within_mid, H.total)
H.to_csv(OUT / "S10h.tsv", sep="\t", index=False)

pd.DataFrame({"family": top10, "mean_share": share_all[top10].values}).to_csv(OUT / "top10.tsv", sep="\t", index=False)
print("families:", len(FAMS), "| placeholders:", sum(f.lower().startswith("unclassified") for f in FAMS),
      "| orders:", omap.nunique())
print("opposing maps with positive turnover:", OPP)
