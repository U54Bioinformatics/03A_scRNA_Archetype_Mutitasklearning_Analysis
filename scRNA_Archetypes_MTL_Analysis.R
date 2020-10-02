library(Seurat)
library(ggplot2)

#### MERGE ALL 10X and iCell8 Seurat objects and perform batch-correction with CCA ------------------------
anchors.3 <- FindIntegrationAnchors(object.list = list(seu.obj1, seu.obj2, seu.obj3, 
                                                       seu.obj4, seu.obj5, seu.obj6, 
                                                       seu.obj7, seu.obj8, seu.obj9), dims = 1:5)
integrated.3 <- IntegrateData(anchorset = anchors.3, dims = 1:5)

# Run the standard workflow for visualization and clustering
integrated.3 <- ScaleData(object = integrated.3, verbose = FALSE, assay = "RNA")
integrated.3 <- FindVariableFeatures(object = integrated.3, verbose = FALSE, assay = "RNA")
integrated.3 <- RunPCA(object = integrated.3, npcs = 30, verbose = FALSE, assay = "RNA")
integrated.3 <- RunUMAP(object = integrated.3, reduction = "pca", dims = 1:30, assay = "RNA")
p1 <- DimPlot(object = integrated.3, reduction = "umap", group.by = "Sample")

# Run the standard workflow for visualization and clustering
DefaultAssay(object = integrated.3) <- "integrated"
integrated.3 <- ScaleData(object = integrated.3, verbose = FALSE)
integrated.3 <- RunPCA(object = integrated.3, npcs = 30, verbose = FALSE)
integrated.3 <- RunUMAP(object = integrated.3, reduction = "pca", dims = 1:30)
p2 <- DimPlot(object = integrated.3, reduction = "umap", group.by = "Sample")
p2 <- DimPlot(object = integrated.3, reduction = "umap", group.by = "hpca_main_type.x")

saveRDS(integrated.3, file="./Seurat_CCAintegrated_iCell8and10X.RDS") # Save the merged matrix of counts

# Extract the first 5 principal components from the Seurat object for archetype analysis 
integrated.10pcs <- FetchData(integrated.3, vars=c("PC_1", "PC_2", "PC_3", "PC_4", "PC_5", 
                                                   "PC_6", "PC_7", "PC_8", "PC_9", "PC_10"))

#### Determine the number of archetypes required to fit a polytope that encloses the gene expression data 
#### projected on the first 2 PCs -----------------------------------------------------------------------

library(ParetoTI) 
library(ggplot2)
library(RColorBrewer)

# How many archetypes are needed to enclose the data? 
arc_ks_t = k_fit_pch(data = t(integrated.10pcs), 
                     ks = 3:8, check_installed = TRUE, bootstrap = FALSE,
                     bootstrap_N = 10, sample_prop = 0.65, 
                     bootstrap_type = "s", seed = 345, simplex = FALSE, var_in_dims = FALSE,
                     normalise_var = TRUE)

p1 <- plot_arc_var(arc_ks_t, type = "varexpl", point_size = 2, line_size = 1.5) + ylim(0, 1) + theme_classic(base_size = 8)
p2 <- plot_arc_var(arc_ks_t, type = "res_varexpl", point_size = 2, line_size = 1.5) + theme_classic(base_size = 8)
p3 <- plot_arc_var(arc_ks_t, type = "t_ratio", point_size = 2, line_size = 1.5) + theme_classic(base_size = 8)

pdf("Integrated_Arcs_Preliminary.pdf", width=3, height=6)
gridExtra::grid.arrange(p1, p2, p3)
dev.off()

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

#### Selecting phenotypes enriched at the archetypes 
#### Multi-task learning with multiple gaussian models with group-lasso penalty ----------------------------------
library(glmnet)
# select hallmark pathway ssGSEA scores for analysis
Hpath <- grep("HALLMARK", colnames(ss.merged)) ## ss.merged is the matrix of ssGSEA scores for all cells

X <- as.matrix(ss.merged[, ..Hpath])
rownames(X) <- ss.merged$Cell.ID
X <- X[intersect(ss.merged$Cell.ID, colnames(arc_data_t$S)), ]

Y <- t(arc_data_t$S[, intersect(ss.merged$Cell.ID, colnames(arc_data_t$S))])
cv.fit1 = cv.glmnet(X, Y, family = "mgaussian", alpha=1)

l1se <- match(cv.fit1$lambda.1se, cv.fit1$lambda)
write.csv(data.frame("A1"=cv.fit1$glmnet.fit$beta$y1[, l1se], 
                     "A2"=cv.fit1$glmnet.fit$beta$y2[, l1se],
                     "A3"=cv.fit1$glmnet.fit$beta$y3[, l1se]),
          file="./../Archetype_Seurat/Integrated_10PCs_mgauss_coefs.csv")

integrated_mgauss_coefs <- data.frame("A1"=cv.fit1$glmnet.fit$beta$y1[, l1se], 
                                      "A2"=cv.fit1$glmnet.fit$beta$y2[, l1se],
                                      "A3"=cv.fit1$glmnet.fit$beta$y3[, l1se])

hall.fit <- glmnet(X, Y, family = "mgaussian", alpha=1)
saveRDS(hall.fit, file="Hallmark_MTL_Fit.RDS")
saveRDS(Y, file="~/Desktop/Ovarian/Ovarian/Archetype_Seurat/Old_Integrated_Arc_Scores.RDS")

