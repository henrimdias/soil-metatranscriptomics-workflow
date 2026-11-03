"""Shared paths, loaders and the GSEA engine.

GSEA engine:
weighted (p = 1) KS running sum on the DESeq2 Wald statistic, gene-label
permutation, NES = ES / mean(|null ES| of the same sign), one-sided p in the
direction of ES, leading edge up to the running-sum peak, KOs with identical
18-sample count vectors collapsed to one feature whose set membership is the
union of its members'. Only the gene-set collection changes between analyses.
"""
import os  # paths are read from environment variables (see README)
import pathlib
from collections import defaultdict

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

UP = pathlib.Path(os.path.join(os.environ.get('MTX_ROOT', '.'), ''))
V2 = UP / "results"
IN = V2 / "ko_matrices"
DE = V2 / "deseq2"
SCR = pathlib.Path(os.environ.get('MTX_REF', os.path.join(os.environ.get('MTX_ROOT', '.'), 'reference')))
KO14 = SCR / "kegg" / "KO1-4.txt"
MODTAB = SCR / "kpc" / "ex" / "kegg_pathways_completeness" / "pathways_data" / "modules_table.tsv"
OUT = pathlib.Path(os.environ.get('MTX_OUT', 'out')); OUT.mkdir(parents=True, exist_ok=True)

GROUPS = ["C24", "G24", "C48", "G48", "C72", "G72"]
SAMPLES = [f"{t}_{h}h_{r}" for t in "CG" for h in (24, 48, 72) for r in (1, 2, 3)]


def group_of(sample):
    t, h, _ = sample.split("_")
    return f"{t}{h[:-1]}"


# coefficient order [Intercept, trt_glucose, tp48, tp72, trt:tp48, trt:tp72]
CONTRASTS = {
    "glucose_G24_vs_C24": [0, 1, 0, 0, 0, 0],
    "glucose_G48_vs_C48": [0, 1, 0, 0, 1, 0],
    "glucose_G72_vs_C72": [0, 1, 0, 0, 0, 1],
    "time_C48_vs_C24": [0, 0, 1, 0, 0, 0],
    "time_C72_vs_C24": [0, 0, 0, 1, 0, 0],
    "time_C72_vs_C48": [0, 0, -1, 1, 0, 0],
    "time_G48_vs_G24": [0, 0, 1, 0, 1, 0],
    "time_G72_vs_G24": [0, 0, 0, 1, 0, 1],
    "time_G72_vs_G48": [0, 0, -1, 1, -1, 1],
    "interaction_48_vs_24": [0, 0, 0, 0, 1, 0],
    "interaction_72_vs_24": [0, 0, 0, 0, 0, 1],
    "interaction_72_vs_48": [0, 0, 0, 0, -1, 1],
}
LABEL = {
    "glucose_G24_vs_C24": "G24 vs C24", "glucose_G48_vs_C48": "G48 vs C48",
    "glucose_G72_vs_C72": "G72 vs C72", "time_C48_vs_C24": "C48 vs C24",
    "time_C72_vs_C24": "C72 vs C24", "time_C72_vs_C48": "C72 vs C48",
    "time_G48_vs_G24": "G48 vs G24", "time_G72_vs_G24": "G72 vs G24",
    "time_G72_vs_G48": "G72 vs G48", "interaction_48_vs_24": "(G48-C48)-(G24-C24)",
    "interaction_72_vs_24": "(G72-C72)-(G24-C24)", "interaction_72_vs_48": "(G72-C72)-(G48-C48)",
}


def ko_names():
    n = pd.read_csv(IN / "ko_names_full.tsv", sep="\t").set_index("KO").gene_name
    h = pd.read_csv(KO14, sep="\t").drop_duplicates("KO").set_index("KO").KoDescription
    return n.combine_first(h)


def symbol(name):
    if not isinstance(name, str):
        return ""
    return name.split(";")[0].split(",")[-1].strip()


def equivalence_groups(counts):
    """Collapse KOs with identical count vectors (rows of `counts`)."""
    key = counts.apply(lambda r: tuple(r), axis=1)
    members = defaultdict(list)
    for ko, k in key.items():
        members[k].append(ko)
    ko2rep, rep2mem = {}, {}
    for mem in members.values():
        mem = sorted(mem)
        rep2mem[mem[0]] = mem
        for ko in mem:
            ko2rep[ko] = mem[0]
    return ko2rep, rep2mem


def enrichment_score(pos, w_hits, n, k, nr):
    cum = np.cumsum(w_hits) / nr
    miss = (pos - np.arange(k)) / (n - k)
    after = cum - miss
    before = np.concatenate(([0.0], cum[:-1])) - miss
    hi, lo = after.max(), before.min()
    return hi if hi > -lo else lo


def leading_edge(pos, w_hits, n, k, es, genes):
    cum = np.cumsum(w_hits) / w_hits.sum()
    miss = (pos - np.arange(k)) / (n - k)
    if es >= 0:
        peak = int(np.argmax(cum - miss))
        return genes[pos[:peak + 1]]
    peak = int(np.argmin(np.concatenate(([0.0], cum[:-1])) - miss))
    return genes[pos[peak:]]


def gsea_one(stat, sets_rep, nperm, rng, min_set=10, max_set=500):
    """stat: Series indexed by representative KO (already deduplicated).
    sets_rep: dict name -> iterable of representative KOs."""
    st = stat.dropna().sort_values(ascending=False)
    genes = st.index.to_numpy(); w = np.abs(st.to_numpy(float)); N = len(genes)
    pos_of = {g: i for i, g in enumerate(genes)}
    idx = {}
    for name, mem in sets_rep.items():
        p = np.sort(np.array([pos_of[m] for m in mem if m in pos_of], dtype=int))
        if min_set <= len(p) <= max_set:
            idx[name] = p
    if not idx:
        return pd.DataFrame()
    obs = {n_: enrichment_score(p, w[p], N, len(p), w[p].sum()) for n_, p in idx.items()}
    null = {n_: np.empty(nperm) for n_ in idx}
    for j in range(nperm):
        perm = rng.permutation(N)
        for n_, p in idx.items():
            pp = np.sort(perm[p]); wp = w[pp]
            null[n_][j] = enrichment_score(pp, wp, N, len(pp), wp.sum())
    rows = []
    for n_, p in idx.items():
        es = obs[n_]; nl = null[n_]
        same = nl[nl >= 0] if es >= 0 else nl[nl < 0]
        den = np.abs(same).mean() if same.size else np.nan
        nes = es / den if den and np.isfinite(den) and den > 0 else np.nan
        pv = ((np.sum(nl >= es) + 1) / (np.sum(nl >= 0) + 1) if es >= 0
              else (np.sum(nl <= es) + 1) / (np.sum(nl < 0) + 1))
        le = leading_edge(p, w[p], N, len(p), es, genes)
        rows.append(dict(set=n_, size=len(p), ES=es, NES=nes, pval=pv,
                         leading_edge_n=len(le), leading_edge=";".join(map(str, le))))
    r = pd.DataFrame(rows)
    r["padj"] = multipletests(r.pval, method="fdr_bh")[1]
    return r


def sets_to_rep(sets_ko, ko2rep):
    """Map KO-level sets onto deduplicated features (union semantics)."""
    out = {}
    for name, kos in sets_ko.items():
        out[name] = sorted({ko2rep[k] for k in kos if k in ko2rep})
    return out


def run_collection(stats_by_contrast, sets_ko, ko2rep, nperm, seed, min_set=10,
                   max_set=500, tag=""):
    """stats_by_contrast: dict contrast -> Series(stat) indexed by KO (all KOs)."""
    rng = np.random.default_rng(seed)
    sets_rep = sets_to_rep(sets_ko, ko2rep)
    res = []
    for ct, st in stats_by_contrast.items():
        st = st[st.index.map(lambda k: ko2rep.get(k) == k)]
        r = gsea_one(st, sets_rep, nperm, rng, min_set, max_set)
        if r.empty:
            continue
        r.insert(0, "contrast", ct)
        r["contrast_label"] = LABEL.get(ct, ct)
        res.append(r)
        if tag:
            print(f"  [{tag}] {ct:24s} {len(r):3d} sets  {int((r.padj < 0.05).sum()):3d} padj<0.05",
                  flush=True)
    return pd.concat(res, ignore_index=True) if res else pd.DataFrame()
