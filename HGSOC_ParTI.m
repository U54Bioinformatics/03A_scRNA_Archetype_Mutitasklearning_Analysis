% Load the data into Matlab from a comma separated value (CSV) file
% The file is a purely numerical matrix, with patients as rows and genes as
% columns
dat = dlmread('~/Desktop/Ovarian/Ovarian/Archetype_Seurat/CCA_CC_mat_nocol.csv', ',');

% The file is formated as samples (i.e. patients) x genes. 

%% We are now ready to perform Pareto Task Inference.
% We use the Sisal algorithm (1), with up to 8 dimensions. We provide the
% discrete patient attributes, and ask ParTI to preliminary booleanize these
% attributes (0). We also pass continuous patient attributes. We pass a boolean 
% matrix specifiying which genes each continuous feature is baesd on (to be used
% in the leave-one-out procedure). We specify that the enrichment analysis 
% will be performed with a bin size of 5%. Finally, the output of the the 
% analysis will be stored in an Excel spreadsheet, under the name 
% 'CancerRNAseq_enrichmentAnalysis_*.csv'.
[arc, arcFinal, pc] = ParTI_lite(dat, 1, 8, [], ...
    [], 0, [], [], [], 0.05, 'HGSOC_enrichmentAnalysis');


%% Finally, we perform the compete analysis, including randomization
% controls and archetype error estimation.
[arc, arcOrig, pc, errs, pval] = ParTI(dat, 1, 8, [], ...
    [], 0, [], [], [], 0.05, 'HGSOC_enrichmentAnalysis');
