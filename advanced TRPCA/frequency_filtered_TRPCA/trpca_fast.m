function [L,S] = trpca_fast(X)

[n1,n2,n3] = size(X);
n12 = min(n1,n2);
X = reshape(X,n1,n2,n3);
lambda = 1/sqrt(n1*n2);
mu = 1e-4;    % need more iterations but the psnr will higher 
rho = 1.2;
tol = 6e-4;
maxiter = 800;

L = zeros(size(X));
S = zeros(size(X));
Y = zeros(size(X));
Lk = zeros(size(X));
mu = 3e-3;

for iter = 1:maxiter
    L = FDFC(X - S + Y/mu);
    S = prox_l1(X - L + Y/mu,lambda/mu);
    Y = Y + mu * (X - L - S);
    mu = rho*mu; 
    
    err = norm(L(:)- Lk(:))/norm(Lk(:));
    if(iter > 5  && err < tol)
        break;
    else
    	disp(['iter ' num2str(iter) ', mu=' num2str(mu) ', err=' num2str(err)]); 
    end
    Lk = L;
end
