function [ up, down ] = naiveTransferOperators( I, numLevels, rowMajor )
%NAIVETRANSFEROPERATORS Get full-weighting up- and down-sample operators.
% rowMajor - true if image indices

if( ~exist('rowMajor','var') )
    rowMajor=false;
end

up = cell(numLevels,1);
down = cell(numLevels,1);

[rows, cols, ~] = size(I);

 downkern = [1/4,1/2,1/4;
             1/2,  1,1/2;
             1/4,1/2,1/4];
 %downkern = downkern/sum(downkern(:));

for level=numLevels:-1:1
    N = rows*cols;
    if( rowMajor && level==numLevels )
        ndx = reshape(1:N, [cols, rows])';
    else
        ndx = reshape(1:N, [rows, cols]);
    end
    smallndx = reshape(1:(floor(rows/2)*floor(cols/2)), [floor(rows/2),floor(cols/2)]);
    ijk = zeros(9*(floor(rows/2)*floor(cols/2)),3);
    ijkndx = 1;

    % Top-left corner
    ijk(ijkndx,:) = [smallndx(1,1), ndx(1,1), 1];
    ijkndx = ijkndx + 1;
    ijk(ijkndx,:) = [smallndx(1,1), ndx(1,2), 1/2];
    ijkndx = ijkndx + 1;
    ijk(ijkndx,:) = [smallndx(1,1), ndx(2,1), 1/2];
    ijkndx = ijkndx + 1;
    ijk(ijkndx,:) = [smallndx(1,1), ndx(2,2), 1/4];
    ijkndx = ijkndx + 1;
    
    % Top-right corner
    ijk(ijkndx,:) = [smallndx(1,floor(cols/2)), ndx(1,cols-2), 1/2];
    ijkndx = ijkndx + 1;
    ijk(ijkndx,:) = [smallndx(1,floor(cols/2)), ndx(1,cols-1), 1];
    ijkndx = ijkndx + 1;
    ijk(ijkndx,:) = [smallndx(1,floor(cols/2)), ndx(1,cols), 1];
    ijkndx = ijkndx + 1;
    ijk(ijkndx,:) = [smallndx(1,floor(cols/2)), ndx(2,cols-2), 1/4];
    ijkndx = ijkndx + 1;
    ijk(ijkndx,:) = [smallndx(1,floor(cols/2)), ndx(2,cols-1), 1/2];
    ijkndx = ijkndx + 1;
    ijk(ijkndx,:) = [smallndx(1,floor(cols/2)), ndx(2,cols), 1/2];
    ijkndx = ijkndx + 1;
    
    % Bottom-left corner
    ijk(ijkndx,:) = [smallndx(floor(rows/2),1), ndx(rows-2,1), 1/2];
    ijkndx = ijkndx + 1;
    ijk(ijkndx,:) = [smallndx(floor(rows/2),1), ndx(rows-2,2), 1/4];
    ijkndx = ijkndx + 1;
    ijk(ijkndx,:) = [smallndx(floor(rows/2),1), ndx(rows-1,1), 1];
    ijkndx = ijkndx + 1;
    ijk(ijkndx,:) = [smallndx(floor(rows/2),1), ndx(rows-1,2), 1/2];
    ijkndx = ijkndx + 1;
    ijk(ijkndx,:) = [smallndx(floor(rows/2),1), ndx(rows,1), 1];
    ijkndx = ijkndx + 1;
    ijk(ijkndx,:) = [smallndx(floor(rows/2),1), ndx(rows,2), 1/2];
    ijkndx = ijkndx + 1;
    
    % Bottom-right corner
    ijk(ijkndx,:) = [smallndx(floor(rows/2),floor(cols/2)), ndx(rows-2,cols-2), 1/4];
    ijkndx = ijkndx + 1;
    ijk(ijkndx,:) = [smallndx(floor(rows/2),floor(cols/2)), ndx(rows-2,cols-1), 1/2];
    ijkndx = ijkndx + 1;
    ijk(ijkndx,:) = [smallndx(floor(rows/2),floor(cols/2)), ndx(rows-2,cols), 1/2];
    ijkndx = ijkndx + 1;
    ijk(ijkndx,:) = [smallndx(floor(rows/2),floor(cols/2)), ndx(rows-1,cols-2), 1/2];
    ijkndx = ijkndx + 1;
    ijk(ijkndx,:) = [smallndx(floor(rows/2),floor(cols/2)), ndx(rows-1,cols-1), 1];
    ijkndx = ijkndx + 1;
    ijk(ijkndx,:) = [smallndx(floor(rows/2),floor(cols/2)), ndx(rows-1,cols), 1];
    ijkndx = ijkndx + 1;
    ijk(ijkndx,:) = [smallndx(floor(rows/2),floor(cols/2)), ndx(rows,cols-2), 1/2];
    ijkndx = ijkndx + 1;
    ijk(ijkndx,:) = [smallndx(floor(rows/2),floor(cols/2)), ndx(rows,cols-1), 1];
    ijkndx = ijkndx + 1;
    ijk(ijkndx,:) = [smallndx(floor(rows/2),floor(cols/2)), ndx(rows,cols), 1];
    ijkndx = ijkndx + 1;
    
    % Left edge
    for v=3:2:rows-2
        ijk(ijkndx,:) = [smallndx((v+1)/2,1), ndx(v-1,1), 1/2];
        ijkndx = ijkndx + 1;
        ijk(ijkndx,:) = [smallndx((v+1)/2,1), ndx(v-1,2), 1/4];
        ijkndx = ijkndx + 1;
        ijk(ijkndx,:) = [smallndx((v+1)/2,1), ndx(v,1), 1];
        ijkndx = ijkndx + 1;
        ijk(ijkndx,:) = [smallndx((v+1)/2,1), ndx(v,2), 1/2];
        ijkndx = ijkndx + 1;
        ijk(ijkndx,:) = [smallndx((v+1)/2,1), ndx(v+1,1), 1/2];
        ijkndx = ijkndx + 1;
        ijk(ijkndx,:) = [smallndx((v+1)/2,1), ndx(v+1,2), 1/4];
        ijkndx = ijkndx + 1;
    end
    
    % Right edge (really, u=cols-1 line)
    for v=3:2:rows-2
        ijk(ijkndx,:) = [smallndx((v+1)/2,floor(cols/2)), ndx(v-1,cols-2), 1/4];
        ijkndx = ijkndx + 1;
        ijk(ijkndx,:) = [smallndx((v+1)/2,floor(cols/2)), ndx(v-1,cols-1), 1/2];
        ijkndx = ijkndx + 1;
        ijk(ijkndx,:) = [smallndx((v+1)/2,floor(cols/2)), ndx(v-1,cols), 1/2];
        ijkndx = ijkndx + 1;
        ijk(ijkndx,:) = [smallndx((v+1)/2,floor(cols/2)), ndx(v,cols-2), 1/2];
        ijkndx = ijkndx + 1;
        ijk(ijkndx,:) = [smallndx((v+1)/2,floor(cols/2)), ndx(v,cols-1), 1];
        ijkndx = ijkndx + 1;
        ijk(ijkndx,:) = [smallndx((v+1)/2,floor(cols/2)), ndx(v,cols), 1];
        ijkndx = ijkndx + 1;
        ijk(ijkndx,:) = [smallndx((v+1)/2,floor(cols/2)), ndx(v+1,cols-2), 1/4];
        ijkndx = ijkndx + 1;
        ijk(ijkndx,:) = [smallndx((v+1)/2,floor(cols/2)), ndx(v+1,cols-1), 1/2];
        ijkndx = ijkndx + 1;
        ijk(ijkndx,:) = [smallndx((v+1)/2,floor(cols/2)), ndx(v+1,cols), 1/2];
        ijkndx = ijkndx + 1;
    end
    
    % Top edge
    for u=3:2:cols-2
        ijk(ijkndx,:) = [smallndx(1,(u+1)/2), ndx(1,u-1), 1/2];
        ijkndx = ijkndx + 1;
        ijk(ijkndx,:) = [smallndx(1,(u+1)/2), ndx(2,u-1), 1/4];
        ijkndx = ijkndx + 1;
        ijk(ijkndx,:) = [smallndx(1,(u+1)/2), ndx(1,u), 1];
        ijkndx = ijkndx + 1;
        ijk(ijkndx,:) = [smallndx(1,(u+1)/2), ndx(2,u), 1/2];
        ijkndx = ijkndx + 1;
        ijk(ijkndx,:) = [smallndx(1,(u+1)/2), ndx(1,u+1), 1/2];
        ijkndx = ijkndx + 1;
        ijk(ijkndx,:) = [smallndx(1,(u+1)/2), ndx(2,u+1), 1/4];
        ijkndx = ijkndx + 1;
    end
    
    % Bottom edge (really, v=rows-1 line)
    for u=3:2:cols-2
        ijk(ijkndx,:) = [smallndx(floor(rows/2),(u+1)/2), ndx(rows-2,u-1), 1/4];
        ijkndx = ijkndx + 1;
        ijk(ijkndx,:) = [smallndx(floor(rows/2),(u+1)/2), ndx(rows-2,u), 1/2];
        ijkndx = ijkndx + 1;
        ijk(ijkndx,:) = [smallndx(floor(rows/2),(u+1)/2), ndx(rows-2,u+1), 1/4];
        ijkndx = ijkndx + 1;
        ijk(ijkndx,:) = [smallndx(floor(rows/2),(u+1)/2), ndx(rows-1,u-1), 1/2];
        ijkndx = ijkndx + 1;
        ijk(ijkndx,:) = [smallndx(floor(rows/2),(u+1)/2), ndx(rows-1,u), 1];
        ijkndx = ijkndx + 1;
        ijk(ijkndx,:) = [smallndx(floor(rows/2),(u+1)/2), ndx(rows-1,u+1), 1/2];
        ijkndx = ijkndx + 1;
        ijk(ijkndx,:) = [smallndx(floor(rows/2),(u+1)/2), ndx(rows,u-1), 1/2];
        ijkndx = ijkndx + 1;
        ijk(ijkndx,:) = [smallndx(floor(rows/2),(u+1)/2), ndx(rows,u), 1];
        ijkndx = ijkndx + 1;
        ijk(ijkndx,:) = [smallndx(floor(rows/2),(u+1)/2), ndx(rows,u+1), 1/2];
        ijkndx = ijkndx + 1;
    end

    % Interior pixels
    for u=3:2:cols-2
        for v=3:2:rows-2
            winNdx = ndx(v-1:v+1,u-1:u+1);
            i = smallndx((v+1)/2, (u+1)/2);
            
%             tmpkern = double(full(A(ndx(v,u),winNdx(:))) < 0);
%             tmpkern(5) = 1;
%             tmpkern = tmpkern(:) .* downkern(:);
            tmpkern = downkern(:);
            
            ijk(ijkndx:ijkndx+9-1,:) = [repmat([i],[9,1]), winNdx(:), tmpkern(:)];
            ijkndx = ijkndx + 9;
        end
    end
    
    down{level} = sparse(ijk(1:ijkndx-1,1), ijk(1:ijkndx-1,2), ijk(1:ijkndx-1,3), (floor(rows/2)*floor(cols/2)), N);
    up{level} = down{level}';
    
    cols = floor(cols/2);
    rows = floor(rows/2);
end
