function [L,S] = itrpca_bm(X)

[n1, n2, n3] = size(X);
n12 = min(n1,n2);
lambda = 1/sqrt(n3*max(n1,n2));
mu = 1e-4;   
rho = 1.1;
tol = 1e-2;
maxiter = 50;

L = zeros(size(X));
S = zeros(size(X));
Y = zeros(size(X));
Lk = zeros(size(X));
PSNR  = [];
error = [];
mu = 3e-3;

for iter = 1:maxiter    
    L = RSTC(X - S + Y/mu);
	S = prox_l1(X - L + Y/mu,lambda/mu);
    Y = Y + mu * (X - L - S);
    mu = rho*mu; 
    
    err = norm(L(:)- Lk(:))/norm(Lk(:));
    if(iter > maxiter  || err < tol)
        break;
    else
    	disp(['iter ' num2str(iter) ', mu=' num2str(mu) ', err=' num2str(err)]); 
    end
    Lk = L;
end
