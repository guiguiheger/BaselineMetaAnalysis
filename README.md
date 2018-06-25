# BaselineMetaAnalysis
This repository includes scripts and in-house built functions used for meta-analysis of several baseline gene expression data

# BaselineMetaAnalysis
This repository includes scripts and in-house built functions used for meta-analysis of several baseline gene expression data

# The following experiments were used for meta-analysis of baseline matrix generation for OpenTargets

The experiments associated with tissues mainly contributed by GTEx

 * E-MTAB-5214 - GTEx (18K samples)
 * E-MTAB-4344 - ENCODE
 * E-MTAB-3716 - Mammalian Kaessmann 
 * E-MTAB-2836 - Uhlen’s Lab (HPA RNA-seq data) 
 * E-MTAB-513 - Illumina Body Map

 # Blueprint experiments for cell lines
 * E-MTAB-3819
 * E-MTAB-3827 
 * E-MTAB-4754

In this analysis before different baseline expression values were combined across experiments (GTEx and Blueprint), the expression values (raw counts) were aggregated using median statistic within each experiment across technical replicates.

`BaselineMetaAnalysis/OpenTargets/RUV_normalisation_gtex.R`

`BaselineMetaAnalysis/OpenTargets/RUV_normalisation_BluePrint.R`


In order to remove any unwanted variations in the RNA-seq data caused by the batch effects, library preparations or any other nuisance effects, we applied RUVg (Risso, Ngai, Speed, & Dudoit, 2014) using between sample normalisation method. RUVg follows generalized linear model (GLM) that assumes that there are no drastic differences in the gene expression values across samples, and it uses negative control genes (i.e. genes that are non-differentially expressed) as training set to remove any unwanted variations. Coefficient of variation (CV) was estimated for each gene in the dataset across the sample, and 1% quantile threshold was fixed to determine ranked gene list that are least variable. Top 1000 least variable genes were used as training set with default parameters of RUVg to normalise raw counts. It should be noted that no filtration to get rid of lower gene expression values was performed in pre-or-post RUV normalisation step
Mapping of the tissue sample names in the baseline matrix against the tissues in anatomical systems, and all the tissue are correctly mapped were chosen to be final matrix. Samples across the same tissue from different experiments were averaged by median in the final expression matrix, and deciles scores are estimate concomitantly.


The following script is used to combine both GTEx and Blueprints normalized gene expression values
`BaselineMetaAnalysis/OpenTargets/RUV_gtex_blueprint_combine.R`


# References
Risso, D., Ngai, J., Speed, T. P., & Dudoit, S. (2014). Normalization of RNA-seq data using factor analysis of control genes or samples. Nature Biotechnology, 32(9), 896–902. https://doi.org/10.1038/nbt.2931

