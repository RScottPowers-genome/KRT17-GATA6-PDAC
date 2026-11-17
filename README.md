# KRT17/GATA6 tumor-cell categories in PDAC single-cell RNA-seq: gene set enrichment analysis

Code accompanying Delgado-Coka L, et al. "A KRT17/GATA6 hybrid state marks unique vulnerability to MYC" (Scientific Reports, in revision).

## Data

The analysis uses the ALRA-imputed malignant-epithelial single-cell object (77,155 cells) built from Peng et al. 2019 (GSA CRA001160) and Werba et al. 2023 (GEO GSM6204109). Preprocessing of this object (ambient RNA removal, quality control, batch correction, malignant-cell identification and ALRA imputation) is described in Nelson et al., Molecular Oncology 20:1857-1876 (2026), https://doi.org/10.1002/1878-0261.70218. The processed object (NYU_imputed.h5ad, 5.6 GB) is deposited at Zenodo, https://doi.org/10.5281/zenodo.XXXXXXX. The raw data are public at the accessions above.

## Analysis

fourcategory_gsea.py assigns each cell to one of four categories according to whether its imputed KRT17 and GATA6 values are above zero (KRT17+GATA6+, KRT17+GATA6-, KRT17-GATA6+, KRT17-GATA6-). For each category, genes are ranked against all other cells with the Scanpy t-test statistic. Each ranked list is tested with GSEApy preranked GSEA against the MSigDB Hallmark gene sets (1,000 permutations, seed 0, gene set size 15-500).

Run:

    python fourcategory_gsea.py path/to/NYU_imputed.h5ad

Outputs are written to results/.

fourcategory_depmap.py estimates gene dependencies for each of the four categories from DepMap pancreatic cancer cell lines. For each category, mean imputed expression of every gene is computed. Pancreatic cell lines (DepMap lineage "Pancreas") with both expression and CRISPR data are kept, matched between the two DepMap files by depmap_id. Each category's expression profile is correlated with each cell line's profile over shared genes. The five most similar lines are combined, weighted by correlation, into one CRISPR dependency probability per gene. Selectivity is a category's dependency minus the mean of the other three categories.

The DepMap files are not redistributed here. Download from the DepMap portal (https://depmap.org/portal/data_page/?tab=customDownloads), release 26Q1: "Expression (Short-read) Public 26Q1" and "CRISPR Gene Dependency", each with cell-line metadata columns (depmap_id, cell_line_display_name, lineage_1). Run:

    python fourcategory_depmap.py path/to/NYU_imputed.h5ad depmap_expression.csv depmap_crispr_dependency.csv

Outputs, including the list of cell lines used for each category, are written to results_depmap/.

## Note on provenance

The script used to produce the values in the manuscript was not retained. fourcategory_gsea.py reconstructs that analysis on the same input object. results/published_fourcategory_NES.tsv holds the 80 normalized enrichment scores (20 Hallmark sets x 4 categories) reported in the manuscript. The reconstructed analysis reproduces them closely: Pearson r = 0.94 across all 80 values, with the same sign in 95% of them. The script writes this comparison to results/comparison_with_manuscript.tsv.

## Software

Single-cell analyses and figures were generated with Scanpy, including its implementation of PAGA (partition-based graph abstraction; sc.tl.paga) for the trajectory analysis, with GSEApy for gene set enrichment. Parameters are given in the Methods of the manuscript.

## Environment

Python 3.9 or later; package list in environment.yml.

    conda env create -f environment.yml
    conda activate k17g6-gsea

## License

MIT (see LICENSE).

## Contact

R. Scott Powers, Department of Pathology, Stony Brook University.
