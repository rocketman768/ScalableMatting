function [ up, down ] = getTransferOperators( I, labels, numLevels )
%GETTRANSFEROPERATORS Summary of this function goes here
%   Detailed explanation goes here

minInSegment = 6;
up = cell(numLevels,1);
down = cell(numLevels,1);
N = size(I,1)*size(I,2);

I = reshape(I,[N,3]);
labels = reshape(labels, [N,1]);
newLabels = zeros(N,1,'uint32');
ndx = reshape(1:N, size(labels));

minLabel = min(labels);
labels = labels - minLabel + 1;
minLabel = 1;
maxLabel = max(labels);

% COME BACK TO ME!
upijk = zeros(N,3);
downijk = zeros(N,3);

for level=1:numLevels-1
    ijkndx = 1;
    newLabel = 1;
    for label=minLabel:maxLabel
        A = I(labels==label, :);
        segNdx = ndx(labels==label);
        numPix = size(A,1);
        
        % See if pixels are saturated.
        saturated = sum(A>=1,2)>0;
        numSaturated = sum(saturated);
        
        % Cluster is too small, so keep it intact.
        if( numPix < minInSegment )
            upijk(ijkndx,:) = [ newLabel, label, 1];
            downijk(ijkndx,:) = [ label, newLabel, 1];
            ijkndx = ijkndx + 1;
            newLabels(segNdx) = newLabel;
            newLabel = newLabel + 1;
        % Saturation
        elseif( numSaturated >= minInSegment && numSaturated < numPix )
            upijk(ijkndx,:) = [ newLabel, label, 1];
            downijk(ijkndx,:) = [ label, newLabel, 0.5];
            ijkndx = ijkndx + 1;
            newLabels( segNdx(saturated) ) = newLabel;
            newLabel = newLabel + 1;
            
            
            upijk(ijkndx,:) = [ newLabel, label, 1];
            downijk(ijkndx,:) = [ label, newLabel, 0.5];
            ijkndx = ijkndx + 1;
            newLabels( segNdx(~saturated) ) = newLabel;
            newLabel = newLabel + 1;
        % Color line splitting
        else
            mu = sum(A,1)'./numPix;
            [V,D] = eigs(A'*A/numPix - mu*mu',2);
            ratio = D(1,1)/D(2,2);
            
            % Cluster is a blob, so keep it intact.
            %if( ratio < 5 || D(1,1) < 1e-3 )
            if( ratio < 10 || D(1,1) < 1e-3 )
                upijk(ijkndx,:) = [ newLabel, label, 1];
                downijk(ijkndx,:) = [ label, newLabel, 1];
                ijkndx = ijkndx + 1;
                newLabels( segNdx ) = newLabel;
                newLabel = newLabel + 1;
            else
                % Project pixels along dominant color line.
                proj = A*V(:,1) - mu'*V(:,1);
                
                % Split into two halves, positive projection value and
                % negative.
                positive = (proj >= 0);
                
                upijk(ijkndx,:) = [ newLabel, label, 1];
                downijk(ijkndx,:) = [ label, newLabel, 0.5];
                ijkndx = ijkndx + 1;
                newLabels( segNdx(positive) ) = newLabel;
                newLabel = newLabel + 1;
                
                upijk(ijkndx,:) = [ newLabel, label, 1];
                downijk(ijkndx,:) = [ label, newLabel, 0.5];
                ijkndx = ijkndx + 1;
                newLabels( segNdx(~positive) ) = newLabel;
                newLabel = newLabel + 1;
            end
        end
        
    end
    
    % Should it be 'newLabel-1'?
    up{level} = sparse(upijk(1:ijkndx-1,1), upijk(1:ijkndx-1,2), upijk(1:ijkndx-1,3), newLabel-1, maxLabel);
    down{level} = sparse(downijk(1:ijkndx-1,1), downijk(1:ijkndx-1,2), downijk(1:ijkndx-1,3), maxLabel, newLabel-1);
    
    % Update the next level labels.
    labels = newLabels;
    minLabel = 1;
    maxLabel = newLabel-1;
end

% Make the final jump from current number of labels to number of pixels.
ijkndx = 1;
for label=minLabel:maxLabel
    segNdx = ndx(labels==label);
    numPix = size(segNdx,1);
    
    upijk(ijkndx:ijkndx+numPix-1,:) = [ segNdx, repmat(label,[numPix,1]), ones(numPix,1)];
    downijk(ijkndx:ijkndx+numPix-1,:) = [ repmat(label,[numPix,1]), segNdx, ones(numPix,1)/numPix];
    ijkndx = ijkndx + numPix;
end

up{numLevels} = sparse(upijk(1:ijkndx-1,1), upijk(1:ijkndx-1,2), upijk(1:ijkndx-1,3), N, maxLabel);
down{numLevels} = sparse(downijk(1:ijkndx-1,1), downijk(1:ijkndx-1,2), downijk(1:ijkndx-1,3), maxLabel, N);

end

