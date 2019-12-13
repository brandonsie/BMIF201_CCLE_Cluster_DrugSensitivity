
### Identifying patterns in cancer cell line molecular data and implications in drug response prediction
Maha Shady, Greg Brunette, Brandon Sie, Katherine Duchinski

Part of BMIF 201 course project. 

This repo contains resources for a Shiny app to explore drug sensitivity of our identified clusters. 

#### Data:  

* [CCLE](https://portals.broadinstitute.org/ccle) / [DepMap](https://depmap.org/portal/prism/)  
	* RNA-seq   
	* SNV  
	* Mutational signature
	* PRISM drug sensitivity
* [TCGA](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga)  
	* RNA-seq
	* SNV
	* Mutational signature
* [GTEx](https://gtexportal.org/home/datasets)
	* RNA-seq data

#### Clustering: 

#### Drug Sensitivity:

* One-way ANOVA run on PRISM drug screen replicate collapsed log2 fold change vs. DMSO data for each cluster model (4686 drugs, 209 CCLE cell lines (from our 4 tissue types).
* Tukey pairwise comparisons and BH multiple hypothesis correction performed."
* Each ggplot figure represents the drug sensitivity for one drug across the clusters of one model. A figure can be generated if at least one pair of clusters has an adjusted p value < 0.05.

#### Classifier:


#### Links: 

Shiny App: [https://brandonsie.shinyapps.io/BMIF201\_CCLE\_Cluster\_DrugSensitivity](https://brandonsie.shinyapps.io/BMIF201_CCLE_Cluster_DrugSensitivity)  
Github Repo: [https://github.com/brandonsie/BMIF201\_CCLE\_Cluster\_DrugSensitivity](https://github.com/brandonsie/BMIF201_CCLE_Cluster_DrugSensitivity)