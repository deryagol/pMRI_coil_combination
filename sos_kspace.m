function res = sos_kspace(x ,dim, pnorm)
% res = sos(x [,dim, pnorm])
%
% function computes the square root of sum of squares along dimension dim.
% If dim is not specified, it computes it along the last dimension.
%
% (c) Michael Lustig 2009

if nargin < 2
    dim = size(size(x),2);
end

if nargin < 3
    pnorm = 2;
end
x=ifftshift(ifft2(ifftshift(x))); 

res = (sum(abs(x.^pnorm),dim)).^(1/pnorm);
