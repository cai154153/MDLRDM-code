function [W]=MDLRDM(X,lambda1,lambda2,k,alpha)
%When using the code in your research work, please cite 
%our paper
%Fu,Z. et al. (2021) A hierarchical weighted low-rank representation for image clustering and classification. PATTERN RECOGNITION, 112. 
%Jiang,L. et al. (2019) Discovering Cancer Subtypes via an Accurate Fusion Strategy on Multiple Profile Data. FRONTIERS IN GENETICS, 10. 

Z{1} = HWLRR(X{1}, lambda1, lambda2, k, 99);
Z{2} = HWLRR(X{2}, lambda1, lambda2, k, 99);
Z{3} = HWLRR(X{3}, lambda1, lambda2, k, 99);
W = SKF(Z, k, 50,alpha);