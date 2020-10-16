function [L,S] = btrpca(Xn)


% imageName=strcat(num2str(i),'.jpg');
% X = double(imread(imageName));  
% X = X/255;
% maxP = max(abs(X(:)));
% [n1,n2,n3] = size(X);
% Xn = X;
% rhos = 0.1;
% ind = find(rand(n1*n2*n3,1)<rhos);
% Xn(ind) = rand(length(ind),1);
% Xn_size=size(Xn);

% initialize the parameters
k = 0;
K = 1;
rho =1.4;
mu =0.1;
max_mu = 500;
nIter=500;
tol = 1e-2;
psnrMAX=0;
% set block sizes
block_sizes = [24,24,3];

levels = size(block_sizes,1);
ms = block_sizes(:,1);
ns = block_sizes(:,2);
vs = block_sizes(:,3);

N1 = max(ms,ns);
lambda1 = 1;
lambdas = 1/sqrt(N1*vs);
L = zeros(size(Xn)); LU = Xn; 
S = zeros(size(Xn)); Y = zeros(size(Xn));

for it = 1:nIter
    k=k+1;
    Lk = L;
    Sk = S;
    L = blockSVT_tensor(Xn - S -Y/mu, block_sizes,1/ mu);  
    mu = min(rho*mu,max_mu);  
    S = prox_l1(Xn - L-Y/mu,lambdas/mu);
    Y = Y - mu * (Xn - L - S);
    dY = L+S-Xn;
    chgL = max(abs(Lk(:)-L(:)));
    chgS = max(abs(Sk(:)-S(:)));
    chg = max([ chgL chgS max(abs(dY(:))) ]);
     
    if chg < tol
        break;
    end 
end