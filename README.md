# Archetype analysis and multi task learning to determine pathways associated with archetypes for scRNA data

#### Software (tested versions)

System: macOS Catalina, 16GB RAM, 2.8GHz quad-core Intel i7

R (>3.6.3)

R-packages: Seurat (v3.1.5), GSVA (v1.34.0), glmnet (v4.0-2)

[ParetoTI package for R](https://github.com/vitkl/ParetoTI)

Several additional algorithms to find positions of the archetypes are available in the original MATLAB implementation of the Pareto Task Inference method. Also, we use this package to calculate the t-ratios and p-value of the archetypes. 

[ParTI package for MATLAB](https://www.weizmann.ac.il/mcb/UriAlon/download/ParTI)

#### Preparing data for analysis:
The multitask learning analysis performed here requires annotated pre-processed and filtered scRNA-seq data and single-sample gene set enrichment analysis (ssGSEA) scores  corresponding to each cell. 

1. Create Seurat objects for individual scRNA data. Perform QC analyses and filter the data to retain meaningful information, and then export the individual Seurat objects. 

Sample scRNA-seq data and vignette for integration and visualization can be found on the [Seurat website](https://satijalab.org/seurat/v3.2/integration.html).

Also, see this detailed[sample workflow for processing 10X data](https://github.com/U54Bioinformatics/02A_scRNAseq_Seurat) 

2. Calculate ssGSEA scores using the Bioconductor [GSVA package](https://bioconductor.org/packages/release/bioc/html/GSVA.html). 

Alternatively, calculate ssGSEA scores using BETSY with ZINB-WaVE normalized scRNA-seq counts data. [Click here for sample workflow](https://github.com/U54Bioinformatics/02C_scRNAseq_Pathway)


## Follow these steps to determine the number of archetypes and pathway phenotypes associated with archetypes using scRNA-seq data

#### 1. Integrate individual Seurat objects and perform batch correction using CCA. 

```r
# load required packages
library(Seurat)
library(ggplot2)

# Integrated Seurat objects corresponding to individual scRNA-seq datasets   
anchors.3 <- FindIntegrationAnchors(object.list = list(seu.obj1, seu.obj2, seu.obj3), dims = 1:5)
integrated.3 <- IntegrateData(anchorset = anchors.3, dims = 1:5)
```

#### 2. Extract principal components from the batch corrected and merged Seurat object. Try using 5-10 PCs for the analysis.  

```r
# Extract the first 10 principal components from the batch-corrected and integrated Seurat object for archetype analysis 
integrated.10pcs <- FetchData(integrated.3, vars=c("PC_1", "PC_2", "PC_3"))
```

#### 3. Determine the number of archetypes. Fit varying number of polytopes and check proportion of variance explained in the scree plots. Choose "k" or the number of archetypes required to enclose the data based on the elbow of the plot.  

```r
library(ParetoTI) 
library(RColorBrewer)

# How many archetypes are needed to enclose the data? 
arc_ks_t = k_fit_pch(data = t(integrated.10pcs), 
                     ks = 3:8, check_installed = TRUE, bootstrap = FALSE,
                     bootstrap_N = 10, sample_prop = 0.65, 
                     bootstrap_type = "s", seed = 345, simplex = FALSE, var_in_dims = FALSE,
                     normalise_var = TRUE)

# Determine optimal "k" using these plots:
p1 <- plot_arc_var(arc_ks_t, type = "varexpl", point_size = 2, line_size = 1.5) + ylim(0, 1) + theme_classic(base_size = 8)
p2 <- plot_arc_var(arc_ks_t, type = "res_varexpl", point_size = 2, line_size = 1.5) + theme_classic(base_size = 8)
p3 <- plot_arc_var(arc_ks_t, type = "t_ratio", point_size = 2, line_size = 1.5) + theme_classic(base_size = 8)

pdf("Integrated_Arcs_Preliminary.pdf", width=3, height=6)
gridExtra::grid.arrange(p1, p2, p3)
dev.off()

```

The significance and p-value for the selected number of archetypes can be calculated using MATLAB ParTI package.
```matlab
% Load numerical matrix, with patients as rows and genes as
% columns. The matrix was exported from the Seurat scRNA-seq object 
% after CCA normalization and integration 
dat = dlmread('~/CCA_CC_mat_nocol.csv', ',');


%% We are now ready to perform Pareto Task Inference.
% We use the Sisal algorithm (1), with up to 8 dimensions. 
[arc, arcFinal, pc] = ParTI_lite(dat, 1, 8);


%% Finally, we perform the complete analysis, including randomization
% controls and archetype error estimation.
[arc, arcOrig, pc, errs, pval] = ParTI(dat, 1, 8);
```


#### 4: Fit the polytope for the optimum number of archtype. This object will also contain archetype scores for each cell. 

```r
# Fit the polytopes based on optimal number of archetypes and plot 
arc_data_t = fit_pch(t(integrated.10pcs), noc = as.integer(3), delta = 0)

# 2D plots
p1 <- plot_arc(arc_data = arc_data_t, data = t(integrated.10pcs),
               data_lab = apply(arc_data_t$S, 2, max), #adds color to vertex by max archetype score
               data_alpha = 0.75,
               data_size = 2,
               which_dimensions = 1:2) + 
  theme_classic(base_size=18)

library(wesanderson)
COLS <- wes_palette("Darjeeling2", 9, type = c("continuous"))
CO <- rank(unique(integrated.10pcs.annot$New_Patient_ID))
p2 <- plot_arc(arc_data = arc_data_t, data = t(integrated.10pcs),
               data_size = 0.2,
               data_alpha = 1, 
               data_lab = integrated.10pcs.annot$New_Patient_ID,
               colors = c(COLS[CO], "red"),
               which_dimensions = 1:2, 
               text_size = 1) + 
  theme_classic(base_size = 28)


pdf("Integrated_Arcs_Patients.PDF", width=10, height=6)
p2
dev.off()
```

#### 5. Use archetype scores to train a multi-task model (group lasso penalty) with hallmark pahtway enrichment scores or gene expression as predictors. The coefficients for the pathways can be used to identify core phenotypes associated with each archetype. 

```r
library(glmnet)

X <- as.matrix(ssgsea.scores) # ssgsea.scores is matrix of Hallmark ssGSEA pathway enrichment scores calculated using GSVA. Cell IDs are in rownames and pathways in colnames

X <- X[intersect(rownames(ssgsea.scores), colnames(arc_data_t$S)), ]
Y <- t(arc_data_t$S[, intersect(rownames(ssgsea.scores), colnames(arc_data_t$S))])

# Perform cross-validation analysis and fit a group-lasso penalized multitask model

cv.fit1 = cv.glmnet(X, Y, family = "mgaussian", alpha=1)
l1se <- match(cv.fit1$lambda.1se, cv.fit1$lambda) #selects a lambda penalty parameter  

write.csv(data.frame("A1"=cv.fit1$glmnet.fit$beta$y1[, l1se], 
                     "A2"=cv.fit1$glmnet.fit$beta$y2[, l1se],
                     "A3"=cv.fit1$glmnet.fit$beta$y3[, l1se]),
          file="./Integrated_10PCs_mgauss_coefs.csv") # save the coefficients from the cross-validation analysis at the selected lambda

```

