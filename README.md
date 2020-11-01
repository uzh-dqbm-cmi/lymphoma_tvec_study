# lymphoma_tvec_study
Scripts and config files used for the T-VEC lymphoma study

## main.R
This file contains the scripts used to read in and process files found in the `filtered_feature_bc_matrix/` folders (outputs of the 10x `cellranger count` pipeline) and also the `filtered_contig_annotations.csv` files (output by the 10x `cellranger vdj` pipeline). The user has to define the paths to the folders and also to the `cell_type_merging.txt` file which defines the rules for cell type merging.

The script performs the following processes:
- Reads in gene expression matrices and creates `SingleR` and `Seurat` objects out of them.
- Uses the `SingleR` package to perform cell typing and narrows down the range of cell types according to the rules defined in `cell_type_merging.txt`
- Reads in V(D)J sequencing results (contigs as defined by the `cellranger` pipeline) into an R dataframe and calculates basic statistics about immune receptors such as clonotype frequency and makes comparisons over multiple samples possible.
- Uses the `Seurat` package to cluster cells based on gene expression and create plots in order to help identify malignant populations and visualize cell type composition changes.
- Allows for plotting scRNA sequencing gene expression and V(D)J sequencing results together.

## cell_type_merging.txt
This textfile contains the cell type merging rules used in `main.R`

## cell_type_sam.py
This Python3 script processes the output samfile of the `cellranger` pipeline (the actual output is a bamfile which can be converted using the `samtools view` command) and sorts reads based on their barcodes into separate samfiles for each cell type. The script takes csv files containing the read barcodes per cell type.
