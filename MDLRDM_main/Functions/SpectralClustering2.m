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
%�ú�������ͼ���ڽӾ��󣬲�ʹ��Ng��Jordan��Weiss���׾����㷨����ڵ�ľ��ࡣCMat:NxN�ڽӾ���   n�������������    groups:Nά����������ͨ�����׾����õ�N�����N��Ⱥ��������

function groups = SpectralClustering2(CKSym,n)

warning off;
N = size(CKSym,1);
MAXiter = 1000; % Maximum number of iterations for KMeans 
REPlic = 20; % Number of replications for KMeans �����ظ�����

% Normalized spectral clustering according to Ng & Jordan & Weiss
% using Normalized Symmetric Laplacian L = I - D^{-1/2} W D^{-1/2}
% ��һ���Գ�������˹����

DN = diag( 1./sqrt(sum(CKSym)+eps) );%����Ⱦ���
LapN = speye(N) - DN * CKSym * DN;%����������˹����
[uN,sN,vN] = svd(LapN);%��������˹����������ֵ�ֽ� 
kerN = vN(:,N-n+1:N);%ȡ����vN��n����С������ֵ��Ӧ����������������¾���kerN
for i = 1:N
    kerNS(i,:) = kerN(i,:) ./ norm(kerN(i,:)+eps);%�Ծ���kerN��ÿһ������һ��
end
%k-means�㷨��Ŀ������С�����ݵ��������������ĵ��ƽ������֮�͡���ˣ�����Ľ��ȡ���ڳ�ʼ�����ĵ��ѡ��
groups = kmeans(kerNS,n,'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');% groups��N*1���������洢����ÿ����ľ����ţ�singleton����ʾ���ĳ����Ϊ�գ��򽫸ôص����ĵ�����Ϊ��ǰ��Զ�����ݵ�
