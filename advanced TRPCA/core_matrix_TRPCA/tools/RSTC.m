function R = RSTC(X)


[n1,n2,n3] = size(X);
n12 = min(n1,n2);
[U, S, V] = tsvd(X);

for i = 1:n3
    reshapeS(:,i) = diag(S(:,:,i));
end

[u,s,v] = svd(reshapeS);
A = u(:,1) * s(1,1) * v(:,1)';
A = repmat(mean(A(:,1:end),2),1,n3);

for i = 1:n3
    SR(:,:,i) = diag(A(:,i));
end

R = tprod(tprod(U,SR),tran(V));

end