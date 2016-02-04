function [ up, down ] = mergeTransferOperators( I, labels, numLevels )
%MERGETRANSFEROPERATORS Summary of this function goes here
%   Detailed explanation goes here

%minSimilarity = 0.9;
up = cell(numLevels,1);
down = cell(numLevels,1);

[rows, cols] = size(labels);
N = rows*cols;

%I = reshape(I,[N,3]);
labels = reshape(labels, [N,1]);
%newLabels = zeros(N,1,'uint32');
ndx = reshape(1:N, size(labels));

labels = labels - min(labels) + 1;
numLabels = max(labels);

%upijk = zeros(N,3);
downijk = zeros(N,3);

% Make the final jump from current number of labels to number of pixels.
ijkndx = 1;
for label=1:numLabels
    segNdx = ndx(labels==label);
    numPix = size(segNdx,1);
    
    upijk(ijkndx:ijkndx+numPix-1,:) = [ segNdx, repmat(label,[numPix,1]), ones(numPix,1)];
    %downijk(ijkndx:ijkndx+numPix-1,:) = [ repmat(label,[numPix,1]), segNdx, ones(numPix,1)/numPix];
    downijk(ijkndx:ijkndx+numPix-1,:) = [ repmat(label,[numPix,1]), segNdx, ones(numPix,1)];
    ijkndx = ijkndx + numPix;
end

up{numLevels} = sparse(upijk(1:ijkndx-1,1), upijk(1:ijkndx-1,2), upijk(1:ijkndx-1,3), N, numLabels);
down{numLevels} = sparse(downijk(1:ijkndx-1,1), downijk(1:ijkndx-1,2), downijk(1:ijkndx-1,3), numLabels, N);

labels = reshape(labels, [rows,cols]);
for level=numLevels-1:-1:1
    % Merge down 1 level.
    [labels, mapping] = merge(I, labels);
    numLabels = max(max(labels));
    
    up{level} = sparse( [1:length(mapping)]', mapping, ones(length(mapping),1), length(mapping), numLabels );
    
    ijkndx = 1;
    for label = 1:numLabels
        components = find(mapping==label);
        num = length(components);
        %downijk(ijkndx:ijkndx+num-1,:) = [ repmat(label,[num,1]), components, repmat(1/num, [num,1]) ];
        downijk(ijkndx:ijkndx+num-1,:) = [ repmat(label,[num,1]), components, repmat(1, [num,1]) ];
        ijkndx = ijkndx + num;
    end
    
    down{level} = sparse(downijk(1:ijkndx-1,1), downijk(1:ijkndx-1,2), downijk(1:ijkndx-1,3), numLabels, length(mapping));
end

end