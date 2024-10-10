%--------------------------------------------------------------------------
% This function takes an adjacency matrix of a graph and computes the 
% clustering of the nodes using the spectral clustering algorithm of 
% Ng, Jordan and Weiss.
% CMat: NxN adjacency matrix
% n: number of groups for clustering
% groups: N-dimensional vector containing the memberships of the N points 
% to the n groups obtained by spectral clustering
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2012
%--------------------------------------------------------------------------
%该函数采用图的邻接矩阵，并使用Ng、Jordan和Weiss的谱聚类算法计算节点的聚类。CMat:NxN邻接矩阵   n：聚类组的组数    groups:N维向量，包含通过光谱聚类获得的N个点对N个群的隶属度

function groups = SpectralClustering2(CKSym,n)

warning off;
N = size(CKSym,1);
MAXiter = 1000; % Maximum number of iterations for KMeans 
REPlic = 20; % Number of replications for KMeans 聚类重复次数

% Normalized spectral clustering according to Ng & Jordan & Weiss
% using Normalized Symmetric Laplacian L = I - D^{-1/2} W D^{-1/2}
% 归一化对称拉普拉斯矩阵

DN = diag( 1./sqrt(sum(CKSym)+eps) );%计算度矩阵
LapN = speye(N) - DN * CKSym * DN;%计算拉普拉斯矩阵
[uN,sN,vN] = svd(LapN);%对拉普拉斯矩阵做奇异值分解 
kerN = vN(:,N-n+1:N);%取矩阵vN的n个最小的特征值对应的特征向量，组成新矩阵kerN
for i = 1:N
    kerNS(i,:) = kerN(i,:) ./ norm(kerN(i,:)+eps);%对矩阵kerN的每一行做归一化
end
%k-means算法的目标是最小化数据点与其所属簇中心点的平方距离之和。因此，聚类的结果取决于初始簇中心点的选择
groups = kmeans(kerNS,n,'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');% groups是N*1的向量，存储的是每个点的聚类标号；singleton’表示如果某个簇为空，则将该簇的中心点设置为当前最远的数据点
