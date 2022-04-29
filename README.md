# Single-cell Developmental Trajectories Inference
This project analyzed the single-cell developmental trajectories during reprogramming. The scRNA-seq dataset was collected by Schiebinger et al. (2019).


## Introduction
A brief introduction about the folders and files:
* `code/`: implementation code;
    * `proposal.Rmd`: Code for proposal;
    * `model.py`: Code for final slides and final report.
* `data/`: used data downloaded from [Google Drive](https://drive.google.com/file/d/1E494DhIx5RLy0qv_6eWa9426Bfmq28po/view);
    * `anndata_subset.h5ad`: gene expression matrix of randomly selected 10% cells (with metadata);
    * `cell_days.txt`: cell sampling days;
    * `cell_sets.gmt`: major cell sets;
    * `ExprMatrix.h5ad`: full gene expression matrix;
    * `ExprMatrix.var.genes.h5ad`: gene expression matrix only with highly variable genes;
    * `fle_coords.txt`: visualization coordinates;
    * `serum_cell_ids.txt`: cells ID.
* `figures/`: output figures;
    * `proposal/`: figures output from `proposal.Rmd`;
    * `final/`: figures output from `model.py`.
* `final_report.pdf`: final report.
* `final_slides.pdf`: final presentation slides.


## Reproducibility
* `proposal.Rmd` can reproduce all the results in proposal.
* `model.py` can reproduce all the results in final slides and final report. 


## Dependencies
* R
    * Seurat
    * SeuratDisk
    * Matrix
    * dplyr
    * ggplot2
    * ggpubr
* Python
    * numpy
    * pandas
    * matplotlib
    * wot
    * os


