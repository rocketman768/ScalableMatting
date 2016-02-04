function [ alpha ] = directSolve( L, c, b, gamma )
%DIRECTSOLVE Summary of this function goes here
%   Detailed explanation goes here

N = size(L,1);
alpha = (L + gamma*spdiags(c,[0], N, N)) \ (gamma*b);

end

