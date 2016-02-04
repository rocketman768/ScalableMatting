function [ A ] = segmentAdjacency( labels )
%SEGMENTADJACENCY Summary of this function goes here
%   Detailed explanation goes here

[rows,cols] = size(labels);

minLabel = min(min(labels));
labels = labels - minLabel + 1;
numLabels = max(max(labels));

A = zeros(numLabels,numLabels);

for v = 2:rows-1
    for u = 2:cols-1
        if( labels(v,u) ~= labels(v-1,u) )
            A(labels(v,u), labels(v-1,u)) = 1;
            A(labels(v-1,u), labels(v,u)) = 1;
        end
        if( labels(v,u) ~= labels(v+1,u) )
            A(labels(v,u), labels(v+1,u)) = 1;
            A(labels(v+1,u), labels(v,u)) = 1;
        end
        if( labels(v,u) ~= labels(v,u-1) )
            A(labels(v,u), labels(v+1,u)) = 1;
            A(labels(v+1,u), labels(v,u)) = 1;
        end
        if( labels(v,u) ~= labels(v,u+1) )
            A(labels(v,u), labels(v,u+1)) = 1;
            A(labels(v,u+1), labels(v,u)) = 1;
        end
    end
end

end

