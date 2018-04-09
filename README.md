# GSheat
Gsheat is an approach to identify and visualise differential pathway expression of an RNA-seq experiment with several contrasts.

## motivation
To provide a reproducible workflow to generate gsheat figure for the dataset **cybrid.xls.tabular.fmt**, but can be customised to perform different analyses.

## files
The workflow consists of a shell script (**gsheat.sh**) to extract expression values from a gene expression count matrix and generate gene-wise rank files for each sample. Then the script then cycles through each gene set in the GMT file and calculates the rank sum, which is then written to a 3 column table (3col.txt).

Next, the R script (**gsheat.R**) reads the 3col.txt file and performs limma differential analysis. The **contrast.mx** file specifies which contrasts are used in the differential analysis. Next, the top ranked gene sets are extracted based on p-value and plotted as a heatmap.
