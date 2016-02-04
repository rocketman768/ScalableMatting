function [ up, down ] = getTransferOperators_old( I, labels, numLevels )
%GETTRANSFEROPERATORS Summary of this function goes here
%   Detailed explanation goes here

minInSegment = 6;
up = cell(numLevels,1);
down = cell(numLevels,1);
N = size(I,1)*size(I,2);

I = reshape(I,[N,3]);
labels = reshape(labels, [N,1]);
ndx = reshape(1:N, size(labels));

minLabel = min(min(labels));
maxLabel = max(max(labels));
numLabels = maxLabel - minLabel + 1;

% COME BACK TO ME!
upijk = zeros(N,3);
downijk = zeros(N,3);
ijkndx = 1;
cluster = 1;

for label=minLabel:maxLabel
    A = I(labels==label, :);
    segNdx = ndx(labels==label);
    numPix = size(A,1);
    
    % See if pixels are saturated.
    saturated = sum(A>=1,2)>0;
    numSaturated = sum(saturated);
    numUnsaturated = numPix - numSaturated;
    
    % Cluster is too small, so keep it intact.
    if( numPix < minInSegment )
        upijk(ijkndx:ijkndx+numPix-1,:) = [ segNdx, repmat(cluster,[numPix,1]), ones(numPix,1)];
        downijk(ijkndx:ijkndx+numPix-1,:) = [ repmat(cluster,[numPix,1]), segNdx, ones(numPix,1)/numPix];
        ijkndx = ijkndx + numPix;
        cluster = cluster + 1;
    % Saturation
    elseif( numSaturated >= minInSegment )
        upijk(ijkndx:ijkndx+numSaturated-1,:) = [ segNdx(saturated), repmat(cluster,[numSaturated,1]), ones(numSaturated,1)];
        downijk(ijkndx:ijkndx+numSaturated-1,:) = [ repmat(cluster,[numSaturated,1]), segNdx(saturated), ones(numSaturated,1)/numSaturated];
        ijkndx = ijkndx + numSaturated;
        cluster = cluster + 1;
        
        
        upijk(ijkndx:ijkndx+numUnsaturated-1,:) = [ segNdx(~saturated), repmat(cluster,[numUnsaturated,1]), ones(numUnsaturated,1)];
        downijk(ijkndx:ijkndx+numUnsaturated-1,:) = [ repmat(cluster,[numUnsaturated,1]), segNdx(~saturated), ones(numUnsaturated,1)/numUnsaturated];
        cluster = cluster + 1;
        ijkndx = ijkndx + numUnsaturated;
    % Color line splitting
    else
        mu = sum(A,1)'./numPix;
        [V,D] = eigs(A'*A/numPix - mu*mu',2);
        ratio = D(1,1)/D(2,2);
        
        % Cluster is a blob, so keep it intact.
        if( ratio < 5 || D(1,1) < 1e-3 )
            upijk(ijkndx:ijkndx+numPix-1,:) = [ segNdx, repmat(cluster,[numPix,1]), ones(numPix,1)];
            downijk(ijkndx:ijkndx+numPix-1,:) = [ repmat(cluster,[numPix,1]), segNdx, ones(numPix,1)/numPix];
            ijkndx = ijkndx + numPix;
            cluster = cluster + 1;
        else
            % Project pixels along dominant color line.
            proj = A*V(:,1) - mu'*V(:,1);
            
            % Split into two halves, positive projection value and
            % negative.
            positive = (proj >= 0);
            numPositive = sum(positive);
            numNegative = numPix-numPositive;
            
            upijk(ijkndx:ijkndx+numPositive-1,:) = [ segNdx(positive), repmat(cluster,[numPositive,1]), ones(numPositive,1)];
            downijk(ijkndx:ijkndx+numPositive-1,:) = [ repmat(cluster,[numPositive,1]), segNdx(positive), ones(numPositive,1)/numPositive];
            ijkndx = ijkndx + numPositive;
            cluster = cluster + 1;
            
            upijk(ijkndx:ijkndx+numNegative-1,:) = [ segNdx(~positive), repmat(cluster,[numNegative,1]), ones(numNegative,1)];
            downijk(ijkndx:ijkndx+numNegative-1,:) = [ repmat(cluster,[numNegative,1]), segNdx(~positive), ones(numNegative,1)/numNegative];
            ijkndx = ijkndx + numNegative;
            cluster = cluster + 1;
        end
    end
    
end

up{1} = sparse(upijk(1:ijkndx-1,1), upijk(1:ijkndx-1,2), upijk(1:ijkndx-1,3), N, cluster);
down{1} = sparse(downijk(1:ijkndx-1,1), downijk(1:ijkndx-1,2), downijk(1:ijkndx-1,3), cluster, N);

end

