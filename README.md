### Archetype analysis and multi task learning to determine pathways associated with archetypes for scRNA data

#### Required software:
R (>3.6.3)
R-packages: Seurat, glmnet 
[ParetoTI](https://github.com/vitkl/ParetoTI)

#### Follow these steps to determine the number of archetypes and pathway phenotypes associated with archetypes using scRNA-seq data. 

Step 1: Create Seurat objects for individual scRNA data. This has to be done before using the example script. Perform QC analysis and filter the data to retain meaningful information, and then export the individual Seurat objects. [Click here for sample workflow](https://github.com/U54Bioinformatics/02A_scRNAseq_Seurat)

Step 2: Integrate individual Seurat objects and perform batch correction with CCA. 

Step 3: Extract principal components from the batch corrected and merged Seurat object. Try using 5-10 PCs for the analysis.  

Step 4: Determine the number of archetypes. Fit varying number of polytopes and check proportion of variance explained in the scree plots. Choose "k" or the number of archetypes required to enclose the data based on the elbow of the plot.  

Step 5: Fit the polytope for the optimum number of archtype. This will also return archetype scores for each cell. 

Step 6: Use archetype scores to train a multi-task model (group lasso penalty) with hallmark pahtway enrichment scores or gene expression as predictors. The coefficients for the pathways can be used to identify core phenotypes associated with each archetype. 
