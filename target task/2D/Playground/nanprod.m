function [Y] = nanprod(X,dim,varargin)

if nargin < 2
    dimension = 1;
else
    dimension = dim;
end

X(isnan(X)) = 1;

Y = prod(X,dimension);
end