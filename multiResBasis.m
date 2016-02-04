function [ E ] = multiResBasis( labels )
%MULTIRESBASIS Summary of this function goes here
%   Detailed explanation goes here

[h,w] = size(labels);
N = h*w;
minLabel = min(min(labels));
numLabels = max(max(labels)) - minLabel + 1;
ndx = reshape(1:N,[w,h]);

E = sparse(ndx(:),labels(:)-minLabel+1,ones(N,1),N,numLabels);
%E = zeros(N, numLabels);
%E(sub2ind(size(E), ndx(:),labels(:)-minLabel+1)) = 1;

% Normalize E
%E = E/(E'*E);

end

