### Archetype analysis and multi task learning to determine pathways associated with archetypes for scRNA data

#### Software (tested versions):
R (>3.6.3)

R-packages: Seurat (v3.1.5), GSVA (v1.34.0), glmnet (v4.0-2)

[ParetoTI](https://github.com/vitkl/ParetoTI)

#### Preparing data for analysis:
The multitask learning analysis performed here requires annotated pre-processed and filtered scRNA-seq data and single-sample gene set enrichment analysis (ssGSEA) scores  corresponding to each cell. 

Prep-step 1: Create Seurat objects for individual scRNA data. Perform QC analyses and filter the data to retain meaningful information, and then export the individual Seurat objects. [Click here for sample workflow](https://github.com/U54Bioinformatics/02A_scRNAseq_Seurat)

Prep-step 2: Calcular ssGSEA scores using the Bioconductor [GSVA package](https://bioconductor.org/packages/release/bioc/html/GSVA.html). Alternatively, calculate ssGSEA scores using BETSY with ZINB-WaVE normalized scRNA-seq counts data. [Click here for sample workflow](https://github.com/U54Bioinformatics/02C_scRNAseq_Pathway)


#### Follow these steps to determine the number of archetypes and pathway phenotypes associated with archetypes using scRNA-seq data. 

Step 1: Integrate individual Seurat objects and perform batch correction using CCA. 

Step 2: Extract principal components from the batch corrected and merged Seurat object. Try using 5-10 PCs for the analysis.  

Step 3: Determine the number of archetypes. Fit varying number of polytopes and check proportion of variance explained in the scree plots. Choose "k" or the number of archetypes required to enclose the data based on the elbow of the plot.  

Step 4: Fit the polytope for the optimum number of archtype. This will also return archetype scores for each cell. 

Step 5: Use archetype scores to train a multi-task model (group lasso penalty) with hallmark pahtway enrichment scores or gene expression as predictors. The coefficients for the pathways can be used to identify core phenotypes associated with each archetype. 
