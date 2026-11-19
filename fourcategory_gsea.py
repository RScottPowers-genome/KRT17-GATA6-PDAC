# fourcategory_gsea.py
#
# Gene set enrichment analysis of four KRT17/GATA6 tumour-cell categories in human PDAC
# single-cell RNA-seq (Delgado-Coka et al., Scientific Reports).
#
# Input: the ALRA-imputed malignant-epithelial AnnData object described in
# Nelson et al., Mol Oncol 2026 (77,155 cells; Peng CRA001160 + Werba GSM6204109),
# deposited at Zenodo as NYU_imputed.h5ad (https://doi.org/10.5281/zenodo.24322447).
#
# Steps:
#   1. Each cell is assigned to one of four categories by whether its imputed KRT17
#      and GATA6 values are above zero: KRT17+GATA6+, KRT17+GATA6-, KRT17-GATA6+,
#      KRT17-GATA6-.
#   2. For each category, genes are ranked against all other cells with the Scanpy
#      rank_genes_groups Wilcoxon rank-sum statistic (several minutes on ~77k cells).
#   3. Each ranked list is tested with GSEApy preranked GSEA against the MSigDB Hallmark
#      gene sets (Enrichr library MSigDB_Hallmark_2020; downloaded by GSEApy, so an
#      internet connection is needed), 1,000 permutations, seed 0, set size 15-500.
#   4. Results for all four categories are written to one table.
#   5. If results/published_fourcategory_NES.tsv is present, the new NES values are
#      compared with the values reported in the manuscript.
#
# Usage:
#   python fourcategory_gsea.py path/to/imputed_epithelial.h5ad [layer_name]
# layer_name is optional; omit it if the ALRA values are in adata.X.

import sys
import os
import numpy as np
import pandas as pd
import scanpy as sc
import gseapy as gp
from scipy.stats import pearsonr, spearmanr

if len(sys.argv) < 2:
    sys.exit("usage: python fourcategory_gsea.py path/to/imputed_epithelial.h5ad [layer_name]")
H5AD = sys.argv[1]
LAYER = sys.argv[2] if len(sys.argv) > 2 else None

GENE_SETS = "MSigDB_Hallmark_2020"
PERMUTATIONS = 1000
SEED = 0
RANK_METHOD = "wilcoxon"
GROUPS = ["KRT17+GATA6+", "KRT17+GATA6-", "KRT17-GATA6+", "KRT17-GATA6-"]
OUTDIR = "results"
OUT = os.path.join(OUTDIR, "fourcategory_gsea_hallmark.tsv")
PUBLISHED = os.path.join(OUTDIR, "published_fourcategory_NES.tsv")
os.makedirs(OUTDIR, exist_ok=True)

print("scanpy", sc.__version__, "| gseapy", gp.__version__)

# ---------------- load ----------------
adata = sc.read_h5ad(H5AD)
if LAYER is not None:
    adata.X = adata.layers[LAYER]
print(adata)


def gene_values(g):
    v = adata[:, g].X
    v = v.toarray() if hasattr(v, "toarray") else np.asarray(v)
    return v.ravel()


# ---------------- 1. categories ----------------
k17 = gene_values("KRT17") > 0
g6 = gene_values("GATA6") > 0
cat = np.where(k17 & g6, "KRT17+GATA6+",
      np.where(k17 & ~g6, "KRT17+GATA6-",
      np.where(~k17 & g6, "KRT17-GATA6+", "KRT17-GATA6-")))
adata.obs["k17g6"] = pd.Categorical(cat, categories=GROUPS)
print("\nCells per category (total %d):" % adata.n_obs)
print(adata.obs["k17g6"].value_counts())
adata.obs["k17g6"].value_counts().rename_axis("category").rename("n_cells") \
    .to_csv(os.path.join(OUTDIR, "fourcategory_cell_counts.tsv"), sep="\t")

# ---------------- 2. gene ranking ----------------
sc.tl.rank_genes_groups(adata, "k17g6", method=RANK_METHOD, reference="rest",
                        key_added="k17g6_rank")

# ---------------- 3-4. preranked GSEA ----------------
rows = []
for grp in GROUPS:
    rg = sc.get.rank_genes_groups_df(adata, group=grp, key="k17g6_rank")
    rnk = (rg[["names", "scores"]].dropna()
             .drop_duplicates("names")
             .sort_values("scores", ascending=False)
             .reset_index(drop=True))
    res = gp.prerank(rnk=rnk, gene_sets=GENE_SETS, permutation_num=PERMUTATIONS,
                     seed=SEED, min_size=15, max_size=500, threads=4,
                     outdir=None, verbose=False)
    d = res.res2d.copy()
    if "Term" not in d.columns:
        d = d.reset_index().rename(columns={"index": "Term"})
    d["group"] = grp
    rows.append(d[["Term", "ES", "NES", "NOM p-val", "FDR q-val", "group"]])

new = pd.concat(rows, ignore_index=True)
new["NES"] = pd.to_numeric(new["NES"])
new.to_csv(OUT, sep="\t", index=False)
print("\nWrote", OUT)

# ---------------- 5. comparison with manuscript values ----------------
if os.path.exists(PUBLISHED):
    pub = pd.read_csv(PUBLISHED, sep="\t")
    pub["Term"] = pub["Term"].str.strip()
    new["Term"] = new["Term"].str.strip()
    m = pub.merge(new, on=["Term", "group"], suffixes=("_pub", "_new"))
    print("\nMatched %d of %d manuscript rows" % (len(m), len(pub)))
    print("All rows: Pearson r = %.3f, Spearman rho = %.3f, same sign = %.0f%%" % (
        pearsonr(m["NES_pub"], m["NES_new"])[0],
        spearmanr(m["NES_pub"], m["NES_new"])[0],
        100 * np.mean(np.sign(m["NES_pub"]) == np.sign(m["NES_new"]))))
    for grp, s in m.groupby("group"):
        print("  %-14s n=%2d  Pearson r = %.3f  max |dNES| = %.2f  same sign = %.0f%%" % (
            grp, len(s), pearsonr(s["NES_pub"], s["NES_new"])[0],
            np.max(np.abs(s["NES_pub"] - s["NES_new"])),
            100 * np.mean(np.sign(s["NES_pub"]) == np.sign(s["NES_new"]))))
    m[["Term", "group", "NES_pub", "NES_new", "p-val", "NOM p-val"]] \
        .sort_values(["Term", "group"]) \
        .to_csv(os.path.join(OUTDIR, "comparison_with_manuscript.tsv"), sep="\t", index=False)
