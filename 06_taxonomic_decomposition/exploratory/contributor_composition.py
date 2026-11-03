#!/usr/bin/env python3
"""Exploratory — taxon-resolved organisation of transcribed functions:
  Q1  Does the family composition of the transcripts of a given function change
      across the incubation?
  Q2  Does that change differ between glucose and control?
plus a first pass at
  A   decomposition of the change in a function's share of the transcriptome into
      (i) change in family shares of the transcriptome ("contributor turnover") and
      (ii) change in the function's share within each family ("within-family change");
  B   effective number of contributing families, D_eff = exp(H), depth-controlled.

Data: X[family, KO, sample] = KO assignments (not reads; multi-KO reads count once
per KO) from counts_ko_taxonomy_full.tsv, restricted to the families of each domain
and to the KOs of the DESeq2 matrices. Sum over families reproduces
counts_ko_{bacteria,fungi}.tsv exactly (checked, 0 mismatches).

Key control: Bacillaceae falls from ~19 % to ~2 % of bacterial mRNA, so ANY function
Bacillaceae transcribes will show a shift in contributors. Every test is therefore
also run against the background, i.e. after dividing each family's share of the
function by its share of the whole transcriptome (Q1-adj, D_eff ratio).

Run: python3 contributor_composition.py <long.tsv|parquet> <KO1-4.txt> <modules_definition.tsv> <outdir>
"""
import sys, pathlib, itertools
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm

LONG, KO14, MODS, OUT = sys.argv[1:5]
OUT = pathlib.Path(OUT); OUT.mkdir(parents=True, exist_ok=True)
U = pathlib.Path(os.path.join(os.environ.get('MTX_ROOT', '.'), ''))
rng = np.random.default_rng(20261012)
NPERM, NBOOT, NRAR = 9999, 2000, 50

GROUPS = ["C24", "G24", "C48", "G48", "C72", "G72"]
SAMPLES = [f"{t}_{h}h_{r}" for t in "CG" for h in (24, 48, 72) for r in (1, 2, 3)]
grp = pd.Series({s: f"{s[0]}{s.split('_')[1][:-1]}" for s in SAMPLES})
trt = pd.Series({s: s[0] for s in SAMPLES}); tim = pd.Series({s: s.split("_")[1] for s in SAMPLES})

# ------------------------------------------------------------------ function sets
h = pd.read_csv(KO14, sep="\t")
h["pid"] = h.Pathway.str.extract(r"\[PATH:(ko\d+)\]")
maps = h.dropna(subset=["pid"]).groupby("pid").KO.apply(set)
mname = (h.dropna(subset=["pid"]).drop_duplicates("pid").set_index("pid").Pathway
           .str.replace(r" \[PATH:.*", "", regex=True))
mods = pd.read_csv(MODS, sep="\t")
mods = mods[mods.module != "module"]
# pre-declared sets: central C, biosynthesis/growth, N (section 7 of the plan)
PRE = {
    "Glycolysis / Gluconeogenesis": ("central C", maps["ko00010"]),
    "Pentose phosphate pathway": ("central C", maps["ko00030"]),
    "Citrate cycle (TCA cycle)": ("central C", maps["ko00020"]),
    "Pyruvate metabolism": ("central C", maps["ko00620"]),
    "Ribosome": ("biosynthesis", maps["ko03010"]),
    "Val/Leu/Ile biosynthesis": ("biosynthesis", maps["ko00290"]),
    "Arginine biosynthesis": ("biosynthesis", maps["ko00220"]),
    "Peptidoglycan biosynthesis": ("biosynthesis", maps["ko00550"]),
    "Nitrogen metabolism": ("nitrogen", maps["ko00910"]),
    "Ntr regulon (module)": ("nitrogen", set(mods[mods.module.str.startswith("Ntr")].KO)),
}
PRE_ID = {"Glycolysis / Gluconeogenesis": "ko00010", "Pentose phosphate pathway": "ko00030",
          "Citrate cycle (TCA cycle)": "ko00020", "Pyruvate metabolism": "ko00620",
          "Ribosome": "ko03010", "Val/Leu/Ile biosynthesis": "ko00290",
          "Arginine biosynthesis": "ko00220", "Peptidoglycan biosynthesis": "ko00550",
          "Nitrogen metabolism": "ko00910", "Ntr regulon (module)": "module"}

# ------------------------------------------------------------------ cube
L = pd.read_parquet(LONG) if LONG.endswith("parquet") else pd.read_csv(LONG, sep="\t")
fam_dom = {}
for d in ("bacteria", "fungi"):
    for f in pd.read_csv(U / f"inputs/counts_family_{d}.tsv", sep="\t").family:
        fam_dom[f] = d.capitalize()
L["domain"] = L.family.map(fam_dom)
L = L.dropna(subset=["domain"])


def cube(dom):
    ko_univ = pd.read_csv(U / f"results/ko_matrices/counts_ko_{dom.lower()}.tsv",
                          sep="\t", index_col=0).index
    d = L[(L.domain == dom) & L.KO.isin(ko_univ)]
    return d.pivot_table(index=["family", "KO"], columns="sample_id", values="count",
                         aggfunc="sum", fill_value=0)[SAMPLES]


# ------------------------------------------------------------------ stats helpers
Xd = pd.get_dummies(pd.DataFrame({"t": trt, "h": tim}), drop_first=False).astype(float)
def design_terms():
    t = (trt == "G").astype(float).to_numpy()[:, None] - 0.5
    hm = pd.get_dummies(tim).astype(float).to_numpy()
    hm = hm[:, 1:] - hm[:, 1:].mean(0)
    th = t * hm
    one = np.ones((18, 1))
    return one, t, hm, th
ONE, T, H, TH = design_terms()


def hat(M):
    return M @ np.linalg.pinv(M)


# balanced design: sequential SS = marginal SS; project onto each term
P_full = hat(np.hstack([ONE, T, H, TH]))
P_terms = {"treatment": hat(np.hstack([ONE, T])) - hat(ONE),
           "time": hat(np.hstack([ONE, H])) - hat(ONE),
           "interaction": P_full - hat(np.hstack([ONE, T, H]))}
DF = {"treatment": 1, "time": 2, "interaction": 2}
DF_RES = 18 - 6


def permanova(Y, nperm=NPERM):
    """Euclidean PERMANOVA (Y already Hellinger-transformed). Free permutation of samples."""
    Yc = Y - Y.mean(0)
    sst = (Yc ** 2).sum()
    def stats(Yp):
        res = ((Yp - P_full @ Yp) ** 2).sum()
        out = {}
        for k, P in P_terms.items():
            ss = ((P @ Yp) ** 2).sum()
            out[k] = (ss / DF[k]) / (res / DF_RES), ss / sst
        return out
    obs = stats(Yc)
    cnt = {k: 1 for k in P_terms}
    for _ in range(nperm):
        Yp = Yc[rng.permutation(18)]
        s = stats(Yp)
        for k in P_terms:
            cnt[k] += s[k][0] >= obs[k][0]
    return {k: (obs[k][1], obs[k][0], cnt[k] / (nperm + 1)) for k in P_terms}


def top_families(C, frac=0.005):
    """families with >= frac of the domain's KO assignments in at least one sample"""
    tot = C.groupby(level="family").sum()
    share = tot / tot.sum()
    return list(share.index[(share >= frac).any(axis=1)])


def collapse(cf, keep):
    """family x sample counts -> kept families + 'Other'"""
    out = cf.reindex(keep).fillna(0)
    out.loc["Other"] = cf.drop(index=[k for k in keep if k in cf.index]).sum()
    return out


def deff(p):
    p = p[p > 0]
    return float(np.exp(-(p * np.log(p)).sum()))


def rarefied_deff(counts, depth, n=NRAR):
    counts = counts.astype(np.int64)
    vals = []
    for _ in range(n):
        sub = rng.multivariate_hypergeometric(counts, depth)
        vals.append(deff(sub / depth))
    return float(np.mean(vals))


# ------------------------------------------------------------------ main
allrows_q, allrows_dec, allrows_def, ko_rows, fam_within = [], [], [], [], []
import os, pickle
MODE = os.environ.get("MODE", "pre")            # pre | screen | combine
DOMS = os.environ.get("DOMS", "Bacteria,Fungi").split(",")
SL = os.environ.get("SLICE", "0/1"); SI, SN = map(int, SL.split("/"))
for dom in (DOMS if MODE != "combine" else []):
    C = cube(dom)
    fam_tot = C.groupby(level="family").sum()                      # family x sample
    keep = top_families(C)
    W = collapse(fam_tot, keep)                                     # kept + Other
    w = W / W.sum()                                                 # family share of transcriptome
    print(f"{dom}: {C.index.get_level_values(0).nunique()} families, "
          f"{C.index.get_level_values(1).nunique()} KOs, {int(C.values.sum()):,} assignments; "
          f"{len(keep)} families kept (>=0.5% in any sample), Other = "
          f"{(W.loc['Other'] / W.sum()).min():.1%}-{(W.loc['Other'] / W.sum()).max():.1%}", flush=True)
    pd.DataFrame(w).to_csv(OUT / f"family_share_{dom}.tsv", sep="\t")

    # background D_eff of the whole transcriptome (rarefied per function depth below)
    sets = dict(PRE)
    # screen: every official map with >= 10 KOs in this domain's universe
    univ = set(C.index.get_level_values(1))
    for pid, kos in maps.items():
        k = kos & univ
        if len(k) >= 10 and int(pid[2:]) < 4000:
            sets.setdefault(mname[pid], ("screen", kos))

    todo = [k for k, v in sets.items() if (v[0] == "screen") == (MODE == "screen")]
    todo = todo[SI::SN]
    KOLEV = C.index.get_level_values(1)
    for fname in todo:
        cls, kos = sets[fname]
        sub = C[KOLEV.isin(kos)]
        if sub.empty:
            continue
        cf = sub.groupby(level="family").sum()
        depth = cf.sum()
        if depth.min() < 200:                                       # coverage filter
            continue
        F = collapse(cf, keep)
        p = F / F.sum()                                             # contributor composition
        # Q1/Q2 raw
        Y = np.sqrt(p.T.to_numpy())
        pr = permanova(Y, 999 if cls == "screen" else NPERM)
        # Q1/Q2 background-adjusted: specialisation e = p / w, renormalised
        e = (p / w.replace(0, np.nan)).fillna(0)
        e = e / e.sum()
        pa = permanova(np.sqrt(e.T.to_numpy()), 999 if cls == "screen" else NPERM)
        row = dict(domain=dom, function=fname, set_class=cls,
                   pathway=PRE_ID.get(fname, next((i for i, n in mname.items() if n == fname), "")),
                   n_KO=len(kos & univ), min_depth=int(depth.min()), median_depth=int(depth.median()))
        for k in P_terms:
            row[f"raw_R2_{k}"], row[f"raw_F_{k}"], row[f"raw_p_{k}"] = pr[k]
            row[f"adj_R2_{k}"], row[f"adj_F_{k}"], row[f"adj_p_{k}"] = pa[k]
        allrows_q.append(row)

        # ---- A: decomposition S = sum_f w_f * phi_f  (phi_f = share of f's transcripts in F)
        phi = (F / W.replace(0, np.nan)).fillna(0)                 # kept+Other x sample
        S = F.sum() / W.sum()
        if cls != "screen" or True:
            def grp_mean(M, g, idx=None):
                cols = [s for s in SAMPLES if grp[s] == g]
                if idx is not None:
                    cols = [cols[i] for i in idx]
                return M[cols].mean(axis=1)
            PAIRS = [("C24", "C48"), ("C24", "C72"), ("G24", "G48"), ("G24", "G72"),
                     ("C24", "G24"), ("C48", "G48"), ("C72", "G72")]
            nb = NBOOT if cls != "screen" else 200
            for a, b in PAIRS:
                def dec(ia=None, ib=None):
                    wa, wb = grp_mean(w, a, ia), grp_mean(w, b, ib)
                    fa, fb_ = grp_mean(phi, a, ia), grp_mean(phi, b, ib)
                    comp = ((wb - wa) * (fa + fb_) / 2)
                    within = ((wa + wb) / 2 * (fb_ - fa))
                    Sa = (wa * fa).sum()
                    return comp.sum() / Sa, within.sum() / Sa, comp / Sa, within / Sa
                c0, w0, cvec, wvec = dec()
                bs = np.array([dec(rng.integers(0, 3, 3), rng.integers(0, 3, 3))[:2] for _ in range(nb)])
                allrows_dec.append(dict(
                    domain=dom, function=fname, set_class=cls, from_group=a, to_group=b,
                    rel_change_total=c0 + w0, turnover=c0, within=w0,
                    turnover_lo=np.percentile(bs[:, 0], 2.5), turnover_hi=np.percentile(bs[:, 0], 97.5),
                    within_lo=np.percentile(bs[:, 1], 2.5), within_hi=np.percentile(bs[:, 1], 97.5),
                    frac_turnover=abs(c0) / (abs(c0) + abs(w0)) if (abs(c0) + abs(w0)) > 0 else np.nan))
                if cls != "screen":
                    for f in F.index:
                        fam_within.append(dict(domain=dom, function=fname, from_group=a, to_group=b,
                                               family=f, turnover=cvec[f], within=wvec[f]))

        # ---- B: D_eff, rarefied to a common depth, and relative to the background
        dep = int(min(depth.min(), 5000))
        for s in SAMPLES:
            de_f = rarefied_deff(F[s].to_numpy(), dep)
            de_b = rarefied_deff(W[s].to_numpy(), dep)
            allrows_def.append(dict(domain=dom, function=fname, set_class=cls, sample=s,
                                    group=grp[s], treatment=trt[s], time=tim[s], depth=dep,
                                    D_eff=de_f, D_eff_background=de_b, D_eff_ratio=de_f / de_b,
                                    top_family=F[s].idxmax(), top_share=float(F[s].max() / F[s].sum())))

        # ---- per-KO Q1/Q2 for the pre-declared sets
        if cls != "screen":
            for ko in sorted(kos & univ):
                cfk = sub[sub.index.get_level_values(1) == ko].droplevel(1)
                if cfk.empty:
                    continue
                if cfk.sum().min() < 50:
                    continue
                Fk = collapse(cfk, keep); pk = Fk / Fk.sum()
                ek = (pk / w.replace(0, np.nan)).fillna(0); ek = ek / ek.sum()
                r1 = permanova(np.sqrt(pk.T.to_numpy()), 999)
                r2 = permanova(np.sqrt(ek.T.to_numpy()), 999)
                ko_rows.append(dict(domain=dom, function=fname, KO=ko, min_depth=int(cfk.sum().min()),
                                    **{f"raw_R2_{k}": r1[k][0] for k in P_terms},
                                    **{f"raw_p_{k}": r1[k][2] for k in P_terms},
                                    **{f"adj_R2_{k}": r2[k][0] for k in P_terms},
                                    **{f"adj_p_{k}": r2[k][2] for k in P_terms}))
        print(f"  {dom} | {fname[:40]:40s} done", flush=True)

if MODE != "combine":
    with open(OUT / f"part_{MODE}_{'_'.join(DOMS)}_{SI}of{SN}.pkl", "wb") as fh:
        pickle.dump((allrows_q, allrows_dec, allrows_def, ko_rows, fam_within), fh)
    print("partial saved"); sys.exit(0)
for pf in sorted(OUT.glob("part_*.pkl")):
    a, b, c, d, e = pickle.load(open(pf, "rb"))
    allrows_q += a; allrows_dec += b; allrows_def += c; ko_rows += d; fam_within += e
q = pd.DataFrame(allrows_q)
for kind in ("raw", "adj"):
    for k in P_terms:
        q[f"{kind}_q_{k}"] = np.nan
        for (dom, cls), idx in q.groupby(["domain", q.set_class == "screen"]).groups.items():
            q.loc[idx, f"{kind}_q_{k}"] = multipletests(q.loc[idx, f"{kind}_p_{k}"], method="fdr_bh")[1]
q.to_csv(OUT / "Q1Q2_permanova_by_function.tsv", sep="\t", index=False)
pd.DataFrame(allrows_dec).to_csv(OUT / "A_decomposition.tsv", sep="\t", index=False)
pd.DataFrame(fam_within).to_csv(OUT / "A_decomposition_by_family.tsv", sep="\t", index=False)
D = pd.DataFrame(allrows_def)
D.to_csv(OUT / "B_Deff_by_sample.tsv", sep="\t", index=False)
kr = pd.DataFrame(ko_rows)
if len(kr):
    for kind in ("raw", "adj"):
        for k in P_terms:
            kr[f"{kind}_q_{k}"] = kr.groupby("domain")[f"{kind}_p_{k}"].transform(
                lambda s: multipletests(s, method="fdr_bh")[1])
    kr.to_csv(OUT / "Q1Q2_permanova_by_KO.tsv", sep="\t", index=False)

# B: two-way ANOVA on D_eff and on the ratio
brow = []
for (dom, fn, cls), d in D.groupby(["domain", "function", "set_class"]):
    out = dict(domain=dom, function=fn, set_class=cls)
    for y in ("D_eff", "D_eff_ratio"):
        a = anova_lm(smf.ols(f"{y} ~ C(treatment) * C(time)", d).fit(), typ=2)
        for term, lab in [("C(treatment)", "treatment"), ("C(time)", "time"),
                          ("C(treatment):C(time)", "interaction")]:
            out[f"{y}_p_{lab}"] = a.loc[term, "PR(>F)"]
        for g in GROUPS:
            out[f"{y}_{g}"] = d[d.group == g][y].mean()
    brow.append(out)
B = pd.DataFrame(brow)
for y in ("D_eff", "D_eff_ratio"):
    for lab in ("treatment", "time", "interaction"):
        B[f"{y}_q_{lab}"] = np.nan
        for _, idx in B.groupby(["domain", B.set_class == "screen"]).groups.items():
            B.loc[idx, f"{y}_q_{lab}"] = multipletests(B.loc[idx, f"{y}_p_{lab}"], method="fdr_bh")[1]
B.to_csv(OUT / "B_Deff_anova.tsv", sep="\t", index=False)
print("done")
