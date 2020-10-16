function R = FCA(X,k)

X = fft(X,[],3);
n3 = size(X,3);
R = zeros(size(X));

if isscalar(k)
	if k == 1
		R(:,:,k) = X(:,:,1);
    else
		R(:,:,k) = X(:,:,k);
		R(:,:,n3+2-k) = conj(R(:,:,k));
    end
else
	R(:,:,k(1):k(end)) = X(:,:,k(1):k(end));
	R(:,:,n3+2-k(1):-1:n3+2-k(end)) = conj(R(:,:,k(1):k(end)));
end

R = ifft(R,[],3);


















