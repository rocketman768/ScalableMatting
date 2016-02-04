function [ outLabels, mapping ] = merge( I, labels )
%MERGE Merges image segments based on mean color similarity.

minLabel = min(min(labels));
labels = labels - minLabel + 1;
numLabels = max(max(labels));

outLabels = zeros(size(labels));

A = segmentAdjacency(labels);

N = size(I,1)*size(I,2);
I = reshape(I, [N, 3]);
labels = reshape(labels, [N,1]);

B = zeros(numLabels,3);
C = zeros(numLabels,numLabels);
merged = zeros(numLabels,1,'int8');
mapping = zeros(numLabels,1);

for label = 1:numLabels
    D = I(labels==label,:);
    B(label,:) = mean(D,1);
end

for i=1:numLabels
    for j=i+1:numLabels
        if( A(i,j) == 0 )
            continue;
        end
        
        C(i,j) = exp(-norm(B(i,:)-B(j,:)));
    end
end

label = 1;
while( true )
    [val,k] = max(C(:));
    if( val < 0.90 )
        break;
    end
    [i,j] = ind2sub([numLabels,numLabels],k);
    merged(i) = 1;
    merged(j) = 1;
    C(i,:) = 0;
    C(:,i) = 0;
    C(j,:) = 0;
    C(:,j) = 0;
    
    outLabels( labels==i | labels==j ) = label;
    mapping(i) = label;
    mapping(j) = label;
    label = label+1;
end

for i=1:numLabels
    if( merged(i) )
        continue;
    end
    
    outLabels( labels==i ) = label;
    mapping(i) = label;
    label = label+1;
end

end