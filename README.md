MDLRDM is written in the MATLAB programming language. To use, please download the MDLRDM folder and follow the instructions provided in the README.

Files：

MDLRDM.m - The main function.

HWLRR.m - Functions used to calculate similarity matrices.

skf.m - Function for fusing similarity matrices.

bestMap.m - Permute labels of L2 to match L1 as good as possible.

hungarian.m - Solve the Assignment problem using the Hungarian method.

Cal_NMI.m - Program for calculating the Normalized Mutual Information (NMI) between two clusterings.

CalcMetrics.m - Program for calculating the accuracy, Normalized Mutual Information (NMI) and error of clustering results

purity.m - Program for calculating the purity of clustering results.

SpectralClustering2.m - The function of spectral clustering

Estimate_Number_of_Clusters_given_graph.m - Calculate the number of clusters.

fkNN.m - Find the nearst neighbors of each sanples in omics data.

constractmap.m - Returns a matrix that records the column index positions of each row of non-zero elements in the input matrix, represented by setting the values of those positions to 1.

discretisation.m - Convert a continuous Ncut vector into a discrete Ncut vector.

discretisationEigenVectorData.m - Discretizes previously rotated eigenvectors in discretisation

sendknew.m - Hierarchy loops to find the nearest neighbors

EProjSimplex_new.m - A function that solves the minimization problem.


We also provided 4 demos, including 3 non cancer datasets for testing the performance of MDLRDM and 1 cancer datasets：

Demo_3sources.m 

Demo_20newsgroups

Demo_BBC4view_685

Demo_colon




Finally, we provided two benchmark methods for evaluating the results, using the R programming language:

survival.R - Calculate the p-value of the cancer classification results, and plot survival analysis.

richment.R - Evaluated the enrichment of clinical labels across cancer classification results.

To run the benchmark, some configuration is needed, mainly settings paths to the datasets, survival and clinical data. 
The expected directory structure is: for the clinical data, a single directory with all clinical data files. 
The datasets and survival data are organized as follows: a single root directory, with subdirectories for every cancer type. 
The directory for each cancer type contains a file for every omic, and a file called "survival" with survival data. 
Omic files, survival data files and clinical label files should be formatted as the files in http://acgt.cs.tau.ac.il/multi_omic_benchmark/download.html.

All cancer datasets from TCGA can be obtianed in http://acgt.cs.tau.ac.il/multi_omic_benchmark/download.html.
