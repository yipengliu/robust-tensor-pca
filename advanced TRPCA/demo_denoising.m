
%% block TRPCA 
clc,clear all;
addpath('./block_TRPCA');
addpath('./block_TRPCA/tools');
Xn = im2double(imread('noise_image.png'));
[L,E] = btrpca(Xn);
figure,
subplot(1,2,1),imshow(L);
subplot(1,2,2),imshow(Xn);


%% core-matrix TRPCA
clc,clear all;
addpath('./core_matrix_TRPCA');
addpath('./core_matrix_TRPCA/tools');
Xn = im2double(imread('noise_image.png'));
[L,E] = btrpca(Xn);
figure,
subplot(1,2,1),imshow(L);
subplot(1,2,2),imshow(Xn);

%% frequency-filter TRPCA
clc,clear all;
addpath('./frequency_filtered_TRPCA');
Xn = im2double(imread('noise_image.png'));
[n1,n2,n3] = size(Xn);
lambda = 1/sqrt(max(n1,n2)*n3);
opts.tau = [.35 1];
% opts.mu = 5e-3;
% opts.rho = 1.3;
% opts.tol = 1e-3;
% opts.max_iter = 50;
[L,E] = trpca_tnn(Xn,lambda,opts);
figure,
subplot(1,2,1),imshow(L);
subplot(1,2,2),imshow(Xn);



