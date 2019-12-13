#### Instructions

The first plot on the right represents drug sensitivity data (PRISM) for a single drug against ~200 cell lines (CCLE) clustered by one of 7 methods. Four inputs help define what data gets plotted. 

1. **Choose Clustering Model**: Select a model id from 1 to 7, corresponding to one of our GMM models clustering CCLE cell lines by various molecular features. 
	* Model name and description will be displayed below the slider.  
	* The bar plots represent the percentage of samples (CCLE cell lines and TCGA samples) that are assigned to each cluster.
2. **Choose Minimum Significant Pairwise Comparisons**: Filter drug data to only plot drugs for which at least **n** pairwise comparisons between clusters from the chosen model are statistically significant.
3. **Choose MOA**: Select a drug mechanism of action from the drop down menu to narrow down possible drugs to plot.
4. **Choose Drug**: Select a drug from the MOA above for which sufficient significant pairwise comparisons are identified.



