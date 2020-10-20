
%% core-matrix TRPCA 
clc,clear all;
addpath('./core_matrix_TRPCA');
addpath('./core_matrix_TRPCA/tools');
load video.mat;
X = im2double(data);
[n1,n2,n3] = size(X);
[L,E] = itrpca_bm(X);
L = abs(L);
E = abs(E);
implay(L)

%% frequency-filter TRPCA
clc,clear all;
addpath('./frequency_filtered_TRPCA');
load video.mat;
data = im2double(data);
[n1,n2,n3] = size(data);
X = reshape(data,n1,n2,n3);
lambda = 1/sqrt(n3*max(n1,n2));
opts.DUBUG = 1;
opts.tau = zeros(1,ceil((n3+1)/2));
opts.tau(2:end) = inf;
opts.mu = 1e-4;
opts.tol = 1e-6;
opts.rho = 1.2;
opts.max_iter = 500;
opts.DEBUG = 1;
[L,E] = trpca_tnn(X,lambda,opts);
% faster version of frequency-filter TRPCA (without SVDs)
% [L,E] = trpca_fast(X); 
L = abs(L);
E = abs(E);
implay(L)




