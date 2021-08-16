# LUSC_COVID-19_Immune
Identification of immune exhaustion class from LUSC cohort using NMF

#####Step 1): We downloaded LUSC FPKM file from TCGA
gdc-client download -m gdc_manifest_20200604_141046.txt

#####Step 2): We renamed the FPKM file using TCGA sample id and moved the files to one directory.
python rename_fpkm.py gdc_sample_sheet.2020-10-13.tsv 


#####Step 3) : We integrated all FPKM files to one expression matrix.
Rscript TCGA_fpkm_reader.R ./ LUSC

#####Step 4) : Estimating the factorization rank
Rscript NMF.R LUSC_protein_Symbols_TCGA_ALL_Samples.txt 

#####Step 5) : calculate stroma and immune enrichment score and generate heatmap based on 4 clusters
Rscript NMF_score_heatmap.R 4 

#####Step 6) : matched clinical data with NMF clusters.
Rscript nmf_map_clindat.R clinical.tsv NMFconsushctree_cluster.xls

#####Step 7) : Kaplan–Meier survival analysis according to clusters
Rscript survival_all_groups_compare.R nmfclusterClindat1.xls 

#####Step 8) : calculate T cell exhaustion and related signatures enrichment scores, and generate heatmap on the basis of scores
Rscript NMF_score_ssgsea_only_heatmap_for_sigmod_xls.R NMFconsushctree_cluster.xls Exhaustion_TEX_immune_cells_cytokine.xls Ovary_protein_Symbols_TCGA_ALL_Samples.txt



Single cell analysis of COVID-19
##### Single cell clustering  of BALF by Seurat
Rscript SingleCell_UMAP_cluster_step1.R

#####cell type annotation and macrophage re-clustering of single cell analysis  in our study
SingleCell_step2.R 


Generated plots in our study, the codes corresponding to Figures 
Figure_plot.R
