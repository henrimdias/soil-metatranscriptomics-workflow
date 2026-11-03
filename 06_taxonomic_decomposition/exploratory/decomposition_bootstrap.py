#!/usr/bin/env python3
"""Exploratory — decomposition at full family resolution (no 'Other' pool).

S_F(s) = sum_f w_f(s) * phi_fF(s)
  w_f    family share of the domain's KO assignments
  phi_fF share of family f's assignments that fall in function F
Change a -> b between group means, midpoint (Shapley) split, exact:
  dS = sum_f dw_f * mean(phi_f)      ("contributor turnover")
     + sum_f mean(w_f) * dphi_f      ("within-family change")
Both expressed relative to S_F(a). Uncertainty: replicate bootstrap within group
(3 of 3 with replacement, independent per group). Also the difference between the
glucose and control time courses for each component, (G24->G72) - (C24->C72), with
the same bootstrap.

e01 pooled families < 0.5 % into one 'Other' unit; a change inside that pool is
partly turnover among minor families but was counted as 'within'. This script is
the sensitivity check (and the version to report).

Run: python3 decomposition_bootstrap.py <long.parquet> <KO1-4.txt> <modules_definition.tsv> <outdir>
"""
import os  # paths are read from environment variables (see README)
import sys, pathlib
import numpy as np
import pandas as pd

LONG, KO14, MODS, OUT = sys.argv[1:5]
OUT = pathlib.Path(OUT)
U = pathlib.Path(os.path.join(os.environ.get('MTX_ROOT', '.'), ''))
rng = np.random.default_rng(7)
GROUPS = ["C24", "G24", "C48", "G48", "C72", "G72"]
SAMPLES = [f"{t}_{h}h_{r}" for t in "CG" for h in (24, 48, 72) for r in (1, 2, 3)]
grp = {s: f"{s[0]}{s.split('_')[1][:-1]}" for s in SAMPLES}
gcols = {g: np.array([i for i, s in enumerate(SAMPLES) if grp[s] == g]) for g in GROUPS}

h = pd.read_csv(KO14, sep="\t")
h["pid"] = h.Pathway.str.extract(r"\[PATH:(ko\d+)\]")
maps = h.dropna(subset=["pid"]).groupby("pid").KO.apply(set)
mname = (h.dropna(subset=["pid"]).drop_duplicates("pid").set_index("pid").Pathway
           .str.replace(r" \[PATH:.*", "", regex=True))
mods = pd.read_csv(MODS, sep="\t")
PRE = {"ko00010": "Glycolysis / Gluconeogenesis", "ko00030": "Pentose phosphate pathway",
       "ko00020": "Citrate cycle (TCA cycle)", "ko00620": "Pyruvate metabolism",
       "ko03010": "Ribosome", "ko00290": "Val/Leu/Ile biosynthesis",
       "ko00220": "Arginine biosynthesis", "ko00550": "Peptidoglycan biosynthesis",
       "ko00910": "Nitrogen metabolism"}
MODSETS = {m: set(mods[mods.module == m].KO) for m in mods.module.unique() if m != "module"}

L = pd.read_parquet(LONG)
fam_dom = {}
for d in ("bacteria", "fungi"):
    for f in pd.read_csv(U / f"inputs/counts_family_{d}.tsv", sep="\t").family:
        fam_dom[f] = d.capitalize()
L["domain"] = L.family.map(fam_dom)

PAIRS = [("C24", "C48"), ("C24", "C72"), ("G24", "G48"), ("G24", "G72"),
         ("C24", "G24"), ("C48", "G48"), ("C72", "G72")]
rows, famrows = [], []
for dom in ("Bacteria", "Fungi"):
    univ = pd.read_csv(U / f"results/ko_matrices/counts_ko_{dom.lower()}.tsv",
                       sep="\t", index_col=0).index
    d = L[(L.domain == dom) & L.KO.isin(univ)]
    C = d.pivot_table(index=["family", "KO"], columns="sample_id", values="count",
                      aggfunc="sum", fill_value=0)[SAMPLES]
    W = C.groupby(level="family").sum()
    fams = W.index
    w = (W / W.sum()).to_numpy()
    kolev = C.index.get_level_values(1)
    sets = {}
    for pid, kos in maps.items():
        if int(pid[2:]) < 4000 and len(kos & set(univ)) >= 10:
            sets[(pid, mname[pid], "pre" if pid in PRE else "screen")] = kos
    for m, kos in MODSETS.items():
        sets[("module", m, "module")] = kos
    for (pid, name, cls), kos in sets.items():
        F = C[kolev.isin(kos)].groupby(level="family").sum().reindex(fams, fill_value=0)
        if F.sum().min() < 200:
            continue
        phi = np.divide(F.to_numpy(), W.to_numpy(), out=np.zeros(F.shape), where=W.to_numpy() > 0)
        nb = 2000 if cls != "screen" else 500
        boot = {}
        for a, b in PAIRS:
            ia = gcols[a][np.vstack([np.arange(3), rng.integers(0, 3, (nb, 3))])]
            ib = gcols[b][np.vstack([np.arange(3), rng.integers(0, 3, (nb, 3))])]
            wa, wb = w[:, ia].mean(2), w[:, ib].mean(2)
            fa, fb = phi[:, ia].mean(2), phi[:, ib].mean(2)
            Sa = (wa * fa).sum(0)
            tv = (wb - wa) * (fa + fb) / 2 / Sa
            wv = (wa + wb) / 2 * (fb - fa) / Sa
            t, wi = tv.sum(0), wv.sum(0)
            boot[(a, b)] = (t, wi)
            rows.append(dict(domain=dom, pathway=pid, function=name, set_class=cls,
                             from_group=a, to_group=b, S_from=Sa[0], S_to=(wb * fb).sum(0)[0],
                             rel_change_total=t[0] + wi[0],
                             turnover=t[0], turnover_lo=np.percentile(t[1:], 2.5), turnover_hi=np.percentile(t[1:], 97.5),
                             within=wi[0], within_lo=np.percentile(wi[1:], 2.5), within_hi=np.percentile(wi[1:], 97.5),
                             frac_turnover=abs(t[0]) / (abs(t[0]) + abs(wi[0]))))
            if cls != "screen":
                order = np.argsort(-np.abs(tv[:, 0]) - np.abs(wv[:, 0]))[:12]
                for j in order:
                    famrows.append(dict(domain=dom, function=name, from_group=a, to_group=b,
                                        family=fams[j], w_from=wa[j, 0], w_to=wb[j, 0],
                                        phi_from=fa[j, 0], phi_to=fb[j, 0],
                                        turnover=tv[j, 0], within=wv[j, 0]))
        # glucose course minus control course (same bootstrap draws reused per pair)
        for (ga, gb), (ca, cb) in [(("G24", "G48"), ("C24", "C48")), (("G24", "G72"), ("C24", "C72"))]:
            tg, wg = boot[(ga, gb)]; tc, wc = boot[(ca, cb)]
            dt, dw = tg - tc, wg - wc
            rows.append(dict(domain=dom, pathway=pid, function=name, set_class=cls,
                             from_group=f"({ga}->{gb})", to_group=f"-({ca}->{cb})",
                             rel_change_total=dt[0] + dw[0],
                             turnover=dt[0], turnover_lo=np.percentile(dt[1:], 2.5), turnover_hi=np.percentile(dt[1:], 97.5),
                             within=dw[0], within_lo=np.percentile(dw[1:], 2.5), within_hi=np.percentile(dw[1:], 97.5),
                             frac_turnover=np.nan))
    print(dom, "done", flush=True)

R = pd.DataFrame(rows)
R["turnover_ci_excl0"] = (R.turnover_lo > 0) | (R.turnover_hi < 0)
R["within_ci_excl0"] = (R.within_lo > 0) | (R.within_hi < 0)
R.to_csv(OUT / "A_decomposition_fullres.tsv", sep="\t", index=False)
pd.DataFrame(famrows).to_csv(OUT / "A_decomposition_fullres_top_families.tsv", sep="\t", index=False)
print(len(R), "rows")
