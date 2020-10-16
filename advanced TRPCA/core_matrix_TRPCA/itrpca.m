function [L,S] = itrpca(X)

mu = 0.5;
tol = 1e-2;
[n1, n2, n3] = size(X);
n12 = min(n1,n2);
lambda1 = 1/sqrt(max(n12,n3));
lambda2 = 1/sqrt(n3*max(n1,n2));

L = zeros(size(X));
S = zeros(size(X));
Y = zeros(size(X));
err = zeros(size(X));
error = [];
iter = [];
niter = 100;

for it = 1:niter
    iter(end+1) = it;
    L = prox_tnn(X - S + Y/mu, mu, lambda1);
    mu = 0.6*mu; 
    figure(1),imshow3(L(:,:,1:40),[],[4,10]),titlef(it);
    drawnow
    S = prox_l1(X - L + Y/mu, lambda2/mu);
    Y = Y + mu * (X - L - S);
    error = norm(err(:)-L(:))/norm(err(:));
    if(it ~= 1 && norm(L(:)- err(:))/norm(err(:)) < tol)
        break;
    else
        disp(['iter=' num2str(it) ',err=' num2str(norm(err(:)-L(:))/norm(err(:)))])
    end
    err = L;
end
Lhat = L;
Shat = S;