function b = blkproc_tensor(varargin)
%BLKPROC Distinct block processing for image.
%
%   BLKPROC is not recommended.  Use BLOCKPROC instead.
%
%   B = BLKPROC(A,[M N],FUN) processes the image A by applying the function
%   FUN to each distinct M-by-N block of A, padding A with zeros if
%   necessary.  FUN is a FUNCTION_HANDLE that accepts an M-by-N matrix, X,
%   and returns a matrix, vector, or scalar Y:
%
%       Y = FUN(X)
%
%   BLKPROC does not require that Y be the same size as X.  However, B is
%   the same size as A only if Y is the same size as X.
%
%   B = BLKPROC(A,[M N],[MBORDER NBORDER],FUN) defines an overlapping
%   border around the blocks.  BLKPROC extends the original M-by-N blocks by
%   MBORDER on the top and bottom, and NBORDER on the left and right,
%   resulting in blocks of size (M+2*MBORDER)-by-(N+2*NBORDER). BLKPROC pads
%   the border with zeros, if necessary, on the edges of A. FUN should
%   operate on the extended block.
%
%   B = BLKPROC(A,'indexed',...) processes A as an indexed image, padding
%   with zeros if the class of A is logical, uint8 or uint16, or ones if 
%   the class of A is double.
%
%   Class Support
%   -------------
%   The input image A can be of any class supported by FUN. The class of B
%   depends on the class of the output from FUN.
%
%   Examples
%   --------
%   FUN can be a FUNCTION_HANDLE created using @.  This example uses BLKPROC
%   to compute the 2-D DCT of each 8-by-8 block of the input image.
%
%       I = imread('cameraman.tif');
%       fun = @dct2;
%       J = blkproc(I,[8 8],fun);
%       imagesc(J), colormap(hot)
%
%   FUN can also be an anonymous function.  This example uses BLKPROC to set
%   the pixels in each 32-by-32 block to the standard deviation of the
%   elements in that block.
%   
%       I = imread('liftingbody.png');
%       fun = @(x) std2(x)*ones(size(x));
%       I2 = blkproc(I,[32 32],fun);
%       figure, imshow(I), figure, imshow(I2,[])
%
%   See also BESTBLK, BLOCKPROC, COLFILT, FUNCTION_HANDLE, NLFILTER.

%   Copyright 1992-2010 The MathWorks, Inc.

% Obsolete syntax: 
%   B = BLKPROC(A,[M N],FUN,P1,P2,...) passes the additional parameters
%   P1,P2,..., to FUN.
%

%   I/O spec
%   ========
%   A      - can be of any class supported by FUN
%   M, N   - double, integer, positive, real scalars
%   FUN    - FUNCTION_HANDLE (could be created with @)

[a, block, border, fun, params, padval] = parse_inputs(varargin{:});

% Expand A: Add border, pad if size(a) is not divisible by block.
[ma,na,va] = size(a);
mpad = rem(ma,block(1)); if mpad>0, mpad = block(1)-mpad; end%rem(x,y)为求整除x/y的余数
npad = rem(na,block(2)); if npad>0, npad = block(2)-npad; end
vpad = rem(va,block(3)); if vpad>0, vpad = block(3)-vpad; end
if (isa(a, 'uint8'))
    if (padval == 1)
        aa = repmat(uint8(1), ma+mpad+2*border(1),na+npad+2*border(2),va+vpad+2*border(3));
    else
        aa = repmat(uint8(0), ma+mpad+2*border(1),na+npad+2*border(2),va+vpad+2*border(3));
    end
elseif isa(a, 'uint16')
    if (padval == 1)
        aa = repmat(uint16(1), ma+mpad+2*border(1),na+npad+2*border(2),va+vpad+2*border(3));
    else
        aa = repmat(uint16(0), ma+mpad+2*border(1),na+npad+2*border(2),va+vpad+2*border(3));
    end
else
    if (padval == 1)
        aa = ones(ma+mpad+2*border(1),na+npad+2*border(2),va+vpad+2*border(3));
    else
        aa = zeros(ma+mpad+2*border(1),na+npad+2*border(2),va+vpad+2*border(3));
    end
end
aa(border(1)+(1:ma),border(2)+(1:na),border(3)+(1:va)) = a;

%% Process first block.
%
m = block(1) + 2*border(1);
n = block(2) + 2*border(2);
v = block(3) + 2*border(3);
mblocks = (ma+mpad)/block(1);
nblocks = (na+npad)/block(2);
vblocks = (va+vpad)/block(3);
arows = 1:m; acols = 1:n;avols=1:v;
x = aa(arows, acols, avols);
firstBlock = feval(fun,x,params{:});
if (isempty(firstBlock))
  style = 'e'; % empty
  b = [];
elseif (all(size(firstBlock) == size(x)))
  style = 's'; % same
  % Preallocate output.
  if (isa(firstBlock, 'uint8'))
     b = repmat(uint8(0), ma+mpad, na+npad, va+vpad);
  elseif (isa(firstBlock, 'uint16'))
     b = repmat(uint16(0), ma+mpad, na+npad, va+vpad);
  else
     b = zeros(ma+mpad, na+npad, va+vpad);
  end
  brows = 1:block(1);
  bcols = 1:block(2);
  bvols = 1:block(3);
  mb = block(1);
  nb = block(2);
  vb = block(3);
  xrows = brows + border(1);
  xcols = bcols + border(2);
  xvols = bvols + border(3);
  b(brows, bcols, bvols) = firstBlock(xrows, xcols,xvols);
elseif (all(size(firstBlock) == (size(x)-2*border)))
  style = 'b'; % border
  % Preallocate output.
  if (isa(firstBlock, 'uint8'))
      b = repmat(uint8(0), ma+mpad, na+npad, va+vpad);
  elseif (isa(firstBlock, 'uint16'))
      b = repmat(uint16(0), ma+mpad, na+npad, va+vpad);
  else
      b = zeros(ma+mpad, na+npad, va+vpad);
  end
  brows = 1:block(1);
  bcols = 1:block(2);
  bvols = 1:block(3);
  b(brows, bcols, bvols) = firstBlock;
  mb = block(1);
  nb = block(2);
  vb = block(3);
else
  style = 'd'; % different
  [P,Q,V] = size(firstBlock);
  brows = 1:P;
  bcols = 1:Q;
  bvols = 1:V;
  mb = P;
  nb = Q;
  vb = V;
  if (isa(firstBlock, 'uint8'))
      b = repmat(uint8(0), mblocks*P, nblocks*Q, vblocks*V);
  elseif (isa(firstBlock, 'uint16'))
      b = repmat(uint16(0), mblocks*P, nblocks*Q, vblocks*V);
  else
      b = zeros(mblocks*P, nblocks*Q, vblocks*V);
  end
  b(brows, bcols, bvols) = firstBlock;
end
%% wait for modify
[rr,cc,vv] = meshgrid(0:(mblocks-1), 0:(nblocks-1), 0:(vblocks-1));
rr = rr(:);
cc = cc(:);
vv = vv(:);
mma = block(1);
nna = block(2);
vva = block(3);
%%process circularly
for k = 2:length(rr)
  x = aa(rr(k)*mma+arows,cc(k)*nna+acols,vv(k)*vva+avols);
  c = feval(fun,x,params{:});
  if (style == 's')
    b(rr(k)*mb+brows,cc(k)*nb+bcols,vv(k)*vb+bvols) = c(xrows,xcols,xvols);
  elseif (style == 'b')
    b(rr(k)*mb+brows,cc(k)*nb+bcols,vv(k)*vb+bvols) = c;
  elseif (style == 'd')
    b(rr(k)*mb+brows,cc(k)*nb+bcols,vv(k)*vb+bvols) = c;
  end
end

if ((style == 's') || (style == 'b'))
  b = b(1:ma,1:na,1:va);
end

%%%
%%% Function parse_inputs
%%%
function [a, block, border, fun, params, padval] = parse_inputs(varargin)

narginchk(2,Inf);

switch nargin
case 3
    % BLKPROC(A, [m n], 'fun')
    a = varargin{1};
    block = varargin{2};
    border = [0 0 0];
    fun = fcnchk(varargin{3});
    params = cell(0,0,0);
    padval = 0;
    
case 4
    if (strcmp(varargin{2}, 'indexed'))
        % BLKPROC(X, 'indexed', [m n], 'fun')
        a = varargin{1};
        block = varargin{3};
        border = [0 0 0];
        fun = fcnchk(varargin{4});
        params = cell(0,0,0);
        padval = 1;
        
    else
        params = varargin(4);
        [fun,msg] = fcnchk(varargin{3}, length(params));
        if isempty(msg)
            % BLKPROC(A, [m n], 'fun', P1)
            a = varargin{1};
            block = varargin{2};
            border = [0 0 0];
            padval = 0;
            
        else
            % BLKPROC(A, [m n], [mb nb], 'fun')
            a = varargin{1};
            block = varargin{2};
            border = varargin{3};
            fun = fcnchk(varargin{4});
            params = cell(0,0,0);
            padval = 0;
        end
    end
    
otherwise
    if (strcmp(varargin{2}, 'indexed'))
        params = varargin(5:end);
        [fun,msg] = fcnchk(varargin{4},length(params));
        if isempty(msg)
            % BLKPROC(A, 'indexed', [m n], 'fun', P1, ...)
            a = varargin{1};
            block = varargin{3};
            border = [0 0 0];
            padval = 1;
            
        else
            % BLKPROC(A, 'indexed', [m n], [mb nb], 'fun', P1, ...)
            a = varargin{1};
            block = varargin{3};
            border = varargin{4};
            params = varargin(6:end);
            fun = fcnchk(varargin{5},length(params));
            padval = 1;
            
        end
        
    else
        params = varargin(4:end);
        [fun,msg] = fcnchk(varargin{3},length(params));
        if isempty(msg)
            % BLKPROC(A, [m n], 'fun', P1, ...)
            a = varargin{1};
            block = varargin{2};
            border = [0 0 0];
            padval = 0;
            
        else
            % BLKPROC(A, [m n], [mb nb], 'fun', P1, ...)
            a = varargin{1};
            block = varargin{2};
            border = varargin{3};
            params = varargin(5:end);
            fun = fcnchk(varargin{4}, length(params));
            padval = 0;
            
        end
        
    end
end
    
if (islogical(a) || isa(a,'uint8') || isa(a, 'uint16'))
    padval = 0;
end

