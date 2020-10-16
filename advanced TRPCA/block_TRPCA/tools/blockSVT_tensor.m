function Z = blockSVT_tensor( Z, block_size, lambda )
% Block-wise singular value thresholding
% 
% Inputs:
% Z          -    input tensor
% lambda     -    threshold
%


blockSVT_tensorfun = @ (X) SVT_tensor( X, lambda );%定义一个自变量为X的函数，

Z = blkproc_tensor( Z, block_size, blockSVT_tensorfun);
% blkproc not recommended in Matlab, but way faster than blockproc...
