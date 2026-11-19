# fourcategory_depmap.py
#
# Transfer of DepMap CRISPR dependency to the four KRT17/GATA6 tumour-cell categories
# (Delgado-Coka et al., Scientific Reports).
#
# Inputs:
#   1. The ALRA-imputed malignant-epithelial AnnData object (Nelson et al., Mol Oncol 2026;
#      NYU_imputed.h5ad, https://doi.org/10.5281/zenodo.24322447).
#   2. DepMap expression, custom download from the DepMap portal
#      (Expression (Short-read) Public 26Q1), one row per cell line, with columns
#      depmap_id, cell_line_display_name, lineage_1 ... followed by genes.
#   3. DepMap CRISPR gene dependency (probability of dependency, 0-1), custom download from
#      the DepMap portal in the same row-per-cell-line layout.
#
# Steps:
#   1. Cells are assigned to KRT17+GATA6+, KRT17+GATA6-, KRT17-GATA6+, KRT17-GATA6- by
#      whether imputed KRT17 and GATA6 are above zero.
#   2. Mean imputed expression of every gene is computed for each category.
#   3. Pancreatic cell lines (lineage_1 contains "pancreas") with both expression and CRISPR
#      data are kept. Cell lines are matched between the two DepMap files by depmap_id.
#   4. Category profiles are centered and scaled within each category across genes; cell-line
#      profiles are centered and scaled within each gene across lines (as in the original
#      analysis). Each category is correlated with each cell line over shared genes.
#   5. For each category, the five most similar lines (positive correlations only, if at
#      least two are positive) are combined, weighted by correlation, into one dependency
#      probability per gene.
#   6. Selectivity = a category's dependency minus the mean of the other three categories.
#
# Usage:
#   python fourcategory_depmap.py NYU_imputed.h5ad depmap_expression.csv depmap_crispr_dependency.csv

import sys
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
from scipy.stats import zscore

if len(sys.argv) < 4:
    sys.exit("usage: python fourcategory_depmap.py H5AD DEPMAP_EXPRESSION_CSV DEPMAP_CRISPR_CSV")
H5AD, EXPR_CSV, CRISPR_CSV = sys.argv[1], sys.argv[2], sys.argv[3]

GROUPS = ["KRT17+GATA6+", "KRT17+GATA6-", "KRT17-GATA6+", "KRT17-GATA6-"]
TOP_N = 5
OUTDIR = "results_depmap"
os.makedirs(OUTDIR, exist_ok=True)


def to_dense(X):
    return X.toarray() if sp.issparse(X) else np.asarray(X)


def load_depmap(path):
    """Read a DepMap portal custom download; return (genes matrix indexed by depmap_id, metadata)."""
    df = pd.read_csv(path, low_memory=False)
    cols_lower = {c: c.lower() for c in df.columns}
    id_col = next((c for c in df.columns if cols_lower[c] == "depmap_id"), None)
    if id_col is None:
        sys.exit("No depmap_id column in %s" % path)
    meta_cols = [c for c in df.columns
                 if cols_lower[c] in ("depmap_id", "cell_line_display_name")
                 or cols_lower[c].startswith("lineage")]
    df[id_col] = df[id_col].astype(str).str.strip()
    meta = df[meta_cols].set_index(id_col)
    mat = df.drop(columns=meta_cols).set_index(df[id_col])
    mat.index.name = "depmap_id"
    mat.columns = mat.columns.astype(str).str.upper()
    mat = mat.apply(pd.to_numeric, errors="coerce")
    if mat.index.duplicated().any():
        sys.exit("Duplicated depmap_id rows in %s" % path)
    return mat, meta


# ---------------- 1-2. categories and category mean expression ----------------
adata = sc.read_h5ad(H5AD)


def gene_values(g):
    return to_dense(adata[:, g].X).ravel()


k17 = gene_values("KRT17") > 0
g6 = gene_values("GATA6") > 0
cat = np.where(k17 & g6, "KRT17+GATA6+",
      np.where(k17 & ~g6, "KRT17+GATA6-",
      np.where(~k17 & g6, "KRT17-GATA6+", "KRT17-GATA6-")))
adata.obs["k17g6"] = pd.Categorical(cat, categories=GROUPS)
print(adata.obs["k17g6"].value_counts())

group_expr = {}
for g in GROUPS:
    group_expr[g] = to_dense(adata[adata.obs["k17g6"] == g].X).mean(axis=0)
group_expr_df = pd.DataFrame(group_expr, index=adata.var_names.astype(str).str.upper())
group_expr_df = group_expr_df[~group_expr_df.index.duplicated()]

# ---------------- 3. DepMap pancreatic lines, matched by depmap_id ----------------
expr, expr_meta = load_depmap(EXPR_CSV)
crispr, crispr_meta = load_depmap(CRISPR_CSV)

lineage_col = next((c for c in expr_meta.columns if c.lower() == "lineage_1"), None)
if lineage_col is None:
    sys.exit("No lineage_1 column in the expression file")
panc_ids = expr_meta.index[expr_meta[lineage_col].astype(str).str.contains("pancreas", case=False)]
shared_ids = pd.Index(panc_ids).intersection(crispr.index)

print("Pancreatic lines with expression:", len(panc_ids))
print("Lines in CRISPR file:", crispr.shape[0])
print("Pancreatic lines with both (matched by depmap_id):", len(shared_ids))

expr_shared = expr.loc[shared_ids]
crispr_shared = crispr.loc[shared_ids]
names = expr_meta.loc[shared_ids, [c for c in expr_meta.columns
                                   if c.lower() == "cell_line_display_name"]]
names.to_csv(os.path.join(OUTDIR, "cell_lines_used.csv"))

# ---------------- 4. similarity ----------------
shared_genes = sorted(set(group_expr_df.index) & set(expr_shared.columns))
print("Shared genes:", len(shared_genes))

tumor_mat = group_expr_df.loc[shared_genes, GROUPS]          # genes x categories
depmap_mat = expr_shared[shared_genes].T                       # genes x lines

tumor_z = pd.DataFrame(zscore(tumor_mat.values, axis=0, nan_policy="omit"),
                       index=tumor_mat.index, columns=tumor_mat.columns).fillna(0)
depmap_z = pd.DataFrame(zscore(depmap_mat.values, axis=1, nan_policy="omit"),
                        index=depmap_mat.index, columns=depmap_mat.columns).fillna(0)

sim = pd.DataFrame(index=GROUPS, columns=depmap_z.columns, dtype=float)
for g in GROUPS:
    gvec = tumor_z[g].values.astype(float)
    for cl in depmap_z.columns:
        sim.loc[g, cl] = np.corrcoef(gvec, depmap_z[cl].values.astype(float))[0, 1]
sim.to_csv(os.path.join(OUTDIR, "category_to_cellline_similarity.csv"))

# ---------------- 5. dependency transfer ----------------
group_dep = {}
lines_used = []
for g in GROUPS:
    top = sim.loc[g].sort_values(ascending=False).head(TOP_N)
    top_pos = top[top > 0]
    if len(top_pos) >= 2:
        top = top_pos
    if len(top) < 2:
        print("Skipping %s: too few matched lines" % g)
        continue
    for cl, r in top.items():
        lines_used.append({"category": g, "depmap_id": cl, "similarity": r})
    dep_sub = crispr_shared.loc[top.index].dropna(axis=1, how="any")
    w = top.values.astype(float)
    w = w / w.sum()
    group_dep[g] = pd.Series(np.average(dep_sub.to_numpy(dtype=float), axis=0, weights=w),
                             index=dep_sub.columns)
    print("%s: lines used %s" % (g, list(top.index)))

pd.DataFrame(lines_used).to_csv(os.path.join(OUTDIR, "matched_lines_per_category.csv"), index=False)
group_dep_df = pd.DataFrame(group_dep)
group_dep_df.to_csv(os.path.join(OUTDIR, "category_transferred_dependency_probabilities.csv"))
print("Dependency matrix:", group_dep_df.shape)

# ---------------- 6. per-category top and selective dependencies ----------------
for g in group_dep_df.columns:
    safe = g.replace("+", "pos").replace("-", "neg")
    top = (group_dep_df[[g]].rename(columns={g: "dependency_score"})
           .rename_axis("gene").reset_index()
           .sort_values("dependency_score", ascending=False))
    top.to_csv(os.path.join(OUTDIR, "%s_top_dependencies.csv" % safe), index=False)

    others = [x for x in group_dep_df.columns if x != g]
    sel = pd.DataFrame({"gene": group_dep_df.index,
                        "dependency_score": group_dep_df[g].values,
                        "other_groups_mean": group_dep_df[others].mean(axis=1).values})
    sel["selectivity_delta"] = sel["dependency_score"] - sel["other_groups_mean"]
    sel = sel.sort_values(["selectivity_delta", "dependency_score"], ascending=[False, False])
    sel.to_csv(os.path.join(OUTDIR, "%s_selective_dependencies.csv" % safe), index=False)
    if "MYC" in group_dep_df.index:
        print("%s: MYC dependency %.3f, selectivity %.3f, rank %d of %d" % (
            g, group_dep_df.at["MYC", g],
            float(sel.loc[sel["gene"] == "MYC", "selectivity_delta"].iloc[0]),
            int(np.where(sel["gene"].values == "MYC")[0][0]) + 1, len(sel)))

print("Done. Saved to", OUTDIR)
