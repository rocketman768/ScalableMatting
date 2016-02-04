function [ p, pinv ] = constraintDistancePerm( scribs )
%CONSTRAINTDISTANCEPERM Gives permuation for distance from scribble
%constraint
%    p - permutation
% pinv - inverse permutation

scribs = scribs(:,:,1);

D = bwdist(scribs > 0.9 | scribs < 0.1,'euclidean');

[~,p] = sort(D(:));
[~,pinv] = sort(p);

end

