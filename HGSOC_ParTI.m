% Load the data into Matlab from a comma separated value (CSV) file
% The file is a purely numerical matrix, with patients as rows and genes as
% columns. The matrix was exported from the Seurat scRNA-seq object 
% after CCA normalization and integration 
dat = dlmread('~/CCA_CC_mat_nocol.csv', ',');


%% We are now ready to perform Pareto Task Inference.
% We use the Sisal algorithm (1), with up to 8 dimensions. 
[arc, arcFinal, pc] = ParTI_lite(dat, 1, 8);


%% Finally, we perform the complete analysis, including randomization
% controls and archetype error estimation.
[arc, arcOrig, pc, errs, pval] = ParTI(dat, 1, 8);
