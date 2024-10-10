clear
clc
close all
addpath("Functions\")
addpath("Datasets\")
load 20newsgroups.mat
x1=data{1};
x2=data{2};
x3=data{3};
for l = 1 : size(x1,2)
    x1(:,l) = x1(:,l)/norm(x1(:,l)); 
end
for l = 1 : size(x2,2)
    x2(:,l) = x2(:,l)/norm(x2(:,l));
end
for l = 1 : size(x3,2)
    x3(:,l) = x3(:,l)/norm(x3(:,l));
end
X{1} = x1;
X{2} = x2;
X{3} = x3;
K1 = length(unique(truth));
k=floor(size(x1,2)/10);
lambda1=1;
lambda2=0.001;
alpha  = 0.6;
W = MDLRDM(X,lambda1,lambda2,k,alpha);
group1 = SpectralClustering2(W,K1);
[ac, nmi, cnt] = CalcMetrics(truth, group1);
[Pri] = purity(truth, group1);