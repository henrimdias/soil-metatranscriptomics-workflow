#!/usr/bin/env python3
"""Independent re-computation of every number plotted by plot_decomposition_masking.R.

Written from scratch (does not import or copy decomposition_contrasts.py / machinery_signature.py): starts from the raw long table
(counts_ko_taxonomy_full.tsv), rebuilds the family x KO x flask cube, and recomputes
  - the decomposition of every KEGG map and of the six gene units (panel A),
  - within-family and community-level changes with 95% CI and P (panel B),
  - the per-flask signature index and the pooled time means with 95% CI (panel C),
then compares each value with the files the plotting script reads.
It also checks the frozen gene lists against their SHA-256 and explains the KEGG map
count. Exit status 1 if any check fails.

Run: python3 verify_decomposition_plot.py <counts_ko_taxonomy_full.tsv> <counts_family_bacteria.tsv>
         <counts_ko_bacteria.tsv> <KO1-4.txt> <machinery_signature_definition.tsv> <.sha256> <figdata_dir>
"""
import sys, hashlib, pathlib
import numpy as np, pandas as pd
from scipy import stats

RAW, FAMB, KOB, KO14, SIG, SHA, FD = sys.argv[1:8]
FD = pathlib.Path(FD)
fails = []
def check(name, a, b, tol=1e-6):
    ok = np.allclose(np.asarray(a, float), np.asarray(b, float), atol=tol, equal_nan=True)
    print(f"[{'OK ' if ok else 'FAIL'}] {name}" + ("" if ok else f"  {a} vs {b}"))
    if not ok:
        fails.append(name)

# ---------------------------------------------------------------- inputs
assert hashlib.sha256(open(SIG, "rb").read()).hexdigest() == open(SHA).read().split()[0]
print("[OK ] gene lists match their frozen SHA-256")
sig = pd.read_csv(SIG, sep="\t")
raw = pd.read_csv(RAW, sep="\t")
fam = set(pd.read_csv(FAMB, sep="\t").family)
kob = pd.read_csv(KOB, sep="\t", index_col=0)
raw = raw[raw.family.isin(fam) & raw.KO.isin(kob.index)]
samples = list(kob.columns)
cube = raw.pivot_table(index=["family", "KO"], columns="sample_id", values="count",
                       aggfunc="sum", fill_value=0)[samples]
# the cube summed over families must reproduce the DESeq2 count matrix exactly
check("cube summed over families == counts_ko_bacteria.tsv",
      cube.groupby(level="KO").sum().reindex(kob.index, fill_value=0).values.sum(), kob.values.sum(), 0)
assert (cube.groupby(level="KO").sum().reindex(kob.index, fill_value=0).values == kob.values).all()

group = np.array([f"{s[0]}{s.split('_')[1].rstrip('h')}" for s in samples])
GROUPS = ["C24", "G24", "C48", "G48", "C72", "G72"]
assert all((group == g).sum() == 3 for g in GROUPS)
contrast = {"C24": -.5, "G24": -.5, "C72": .5, "G72": .5}          # 24 -> 72 h, both treatments


def decompose(W, F):
    """W: family x flask totals; F: family x flask counts of the function."""
    w = W / W.sum(axis=0)
    phi = np.where(W > 0, F / np.where(W > 0, W, 1), 0.0)
    gm = lambda M: np.column_stack([M[:, group == g].mean(1) for g in GROUPS]).mean(1)
    wbar, fbar = gm(w), gm(phi)
    S = (wbar * fbar).sum()
    idx = {"within": (wbar[:, None] * phi).sum(0), "turnover": (w * fbar[:, None]).sum(0),
           "total": (w * phi).sum(0)}
    out = {}
    for k, v in idx.items():
        cm = {g: v[group == g].mean() for g in GROUPS}
        rss = sum(((v[group == g] - cm[g]) ** 2).sum() for g in GROUPS)
        s2 = rss / (18 - 6)
        est = sum(c * cm[g] for g, c in contrast.items())
        se = np.sqrt(s2 * sum(c * c for c in contrast.values()) / 3)
        t = stats.t.ppf(0.975, 12)
        out[k] = dict(est=est / S, lo=(est - t * se) / S, hi=(est + t * se) / S,
                      p=2 * stats.t.sf(abs(est / se), 12), index=v / S)
    # exactness of the split, on group means with the same weights
    return out


W_all = cube.groupby(level="family").sum()
families = W_all.index
kolev = cube.index.get_level_values("KO")

def unit_counts(kos, keep):
    return cube[kolev.isin(kos)].groupby(level="family").sum().reindex(keep, fill_value=0)

# ---------------------------------------------------------------- panel A, KEGG maps
h = pd.read_csv(KO14, sep="\t")
h["pid"] = h.Pathway.str.extract(r"\[PATH:(ko\d+)\]")
maps = h.dropna(subset=["pid"]).groupby("pid").KO.apply(set)
universe = set(kob.index)
eligible = []
for pid, kos in maps.items():
    if int(pid[2:]) >= 4000 or len(kos & universe) < 10:
        continue
    F = unit_counts(kos, families)
    if F.sum(axis=0).min() < 200:
        continue
    eligible.append(pid)
PRE = {"ko00010", "ko00030", "ko00020", "ko00620", "ko03010", "ko00290", "ko00220", "ko00550", "ko00910"}
print(f"       KEGG maps eligible: {len(eligible)} = {len(PRE & set(eligible))} pre-declared + "
      f"{len(set(eligible) - PRE)} screening maps")
fA = pd.read_csv(FD / "decomposition_maps.tsv", sep="\t")
check("panel A map count", len(fA), len(eligible), 0)
mine = {}
for pid in eligible:
    r = decompose(W_all.values.astype(float), unit_counts(maps[pid], families).values.astype(float))
    mine[pid] = (r["turnover"]["est"], r["within"]["est"])
fA = fA.set_index("pathway")
check("panel A map turnover (all 106)", [mine[p][0] for p in fA.index], fA.turnover.values)
check("panel A map within (all 106)", [mine[p][1] for p in fA.index], fA.within.values)

# ---------------------------------------------------------------- units (A, B, C)
LABEL = {"Ribosomal proteins": "Ribosomal proteins",
         "Chaperones (GroESL / DnaK-DnaJ-GrpE)": "Chaperones (GroEL/DnaK)",
         "Sec / SRP translocon": "Sec / SRP", "F-type ATP synthase": "F-ATP synthase",
         "Cell wall: peptidoglycan": "Peptidoglycan", "Cell wall: DAP branch": "DAP branch"}
units = {LABEL[u]: list(d.KO) for u, d in sig.groupby("unit")}
units["Combined signature (4 units)"] = list(sig[~sig.unit.str.startswith("Cell wall")].KO)
fB = pd.read_csv(FD / "decomposition_units.tsv", sep="\t")
res = {}
for variant, keep in (("baseline", families), ("drop_Bacillaceae", families.drop("Bacillaceae"))):
    W = W_all.loc[keep].values.astype(float)
    for lab, kos in units.items():
        res[(variant, lab)] = decompose(W, unit_counts(kos, keep).values.astype(float))
hl = pd.read_csv(FD / "decomposition_maps_highlight.tsv", sep="\t")
check("A  highlighted maps are KEGG maps present in the grey set", hl.pathway.isin(fA.index).all(), True, 0)
check("A  highlighted map coordinates", [[mine[p][0], mine[p][1]] for p in hl.pathway],
      hl[["turnover", "within"]].values)
check("A  panel contains no gene-list unit (all points are ko maps)", all(str(p).startswith("ko") for p in fA.index), True, 0)
for _, row in fB.iterrows():
    r = res[(row.variant, row.label)][row.component]
    check(f"B  {row.variant[:8]:8s} {row.label[:28]:28s} {row.component:8s} est/lo/hi/p",
          [r["est"], r["lo"], r["hi"], r["p"]], [row.rel_change, row.lo, row.hi, row.p])

fC = pd.read_csv(FD / "signature_per_sample.tsv", sep="\t").set_index("sample")
v = res[("baseline", "Combined signature (4 units)")]["within"]["index"]
check("C  per-flask signature index (18 flasks)", v, fC.loc[samples, "within_index_rel"].values)
mC = pd.read_csv(FD / "signature_time_means.tsv", sep="\t").set_index("time")
for tt in (24, 48, 72):
    x = v[np.isin(group, [f"C{tt}", f"G{tt}"])]
    t = stats.t.ppf(.975, len(x) - 1); se = x.std(ddof=1) / np.sqrt(len(x))
    check(f"C  pooled mean and 95% CI at {tt} h (n = {len(x)})",
          [x.mean(), x.mean() - t * se, x.mean() + t * se], mC.loc[tt, ["mean", "lo", "hi"]].values)

# BH within variant x component across the six units (the combined signature is a single test)
from statsmodels.stats.multitest import multipletests
for variant in ("baseline", "drop_Bacillaceae"):
    for comp in ("within", "total", "turnover"):
        labs = [l for l in units if not l.startswith("Combined")]
        q = multipletests([res[(variant, l)][comp]["p"] for l in labs], method="fdr_bh")[1]
        sub = fB[(fB.variant == variant) & (fB.component == comp)].set_index("label")
        have = [l for l in labs if l in sub.index]
        if have:
            check(f"B  {variant[:8]} {comp}: BH q over the six units", [q[labs.index(l)] for l in have], sub.loc[have, "q"].values)

# ---------------------------------------------------------------- facts quoted in the legend
r = res[("baseline", "Combined signature (4 units)")]
x24, x72 = v[np.isin(group, ["C24", "G24"])], v[np.isin(group, ["C72", "G72"])]
print(f"       combined signature: within {100 * r['within']['est']:.1f} %, community "
      f"{100 * r['total']['est']:.1f} % (P = {r['total']['p']:.3f}); without Bacillaceae: within "
      f"{100 * res[('drop_Bacillaceae', 'Combined signature (4 units)')]['within']['est']:.1f} %, community "
      f"{100 * res[('drop_Bacillaceae', 'Combined signature (4 units)')]['total']['est']:.1f} % "
      f"(P = {res[('drop_Bacillaceae', 'Combined signature (4 units)')]['total']['p']:.1e})")
print(f"       all six 72 h flasks below all six 24 h flasks: {x72.max() < x24.min()} "
      f"(max 72 h {x72.max():.3f}, min 24 h {x24.min():.3f})")
print(f"       KO counts per unit: " + ", ".join(f"{k} {len(vv)}" for k, vv in units.items()))
print("\nALL CHECKS PASSED" if not fails else f"\n{len(fails)} CHECK(S) FAILED: {fails}")
sys.exit(1 if fails else 0)
