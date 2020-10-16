function R = FDFC(X)


[n1,n2,n3] = size(X);
R = zeros(n1,n2,n3);
X = fft(X,[],3);
R(:,:,1) = X(:,:,1);


R = ifft(R,[],3);