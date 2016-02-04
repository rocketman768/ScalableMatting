function [ up, down ] = adaptiveTransferOperators( I, A, numLevels )
%NAIVETRANSFEROPERATORS Summary of this function goes here
%   Detailed explanation goes here

up = cell(numLevels,1);
down = cell(numLevels,1);

[rows, cols, ~] = size(I);

for level=numLevels:-1:1
    N = rows*cols;
    ndx = reshape(1:N, [rows, cols]);
    smallndx = reshape(1:(floor(rows/2)*floor(cols/2)), [floor(rows/2),floor(cols/2)]);
    ijk = zeros(9*(floor(rows/2)*floor(cols/2)),3);
    ijkndx = 1;

    % Top-left corner
%     winNdx = ndx(1:2,1:2);
%     i = smallndx(1,1);
%     ijk(ijkndx:ijkndx+4-1,:) = [repmat([i],[4,1]), winNdx(:), reshape(downkern(2:3,2:3)/sum(sum(downkern(2:3,2:3))),[4,1])];
%     ijkndx = ijkndx + 4;
    ijk(ijkndx,:) = [smallndx(1,1), ndx(1,1), 1];
    ijkndx = ijkndx + 1;
    
    % Top-right corner
%     winNdx = ndx(1:2,cols-1:cols);
%     i = smallndx(1,cols/2);
%     ijk(ijkndx:ijkndx+4-1,:) = [repmat([i],[4,1]), winNdx(:), reshape(downkern(2:3,1:2)/sum(sum(downkern(2:3,1:2))),[4,1])];
%     ijkndx = ijkndx + 4;
    ijk(ijkndx,:) = [smallndx(1,floor(cols/2)), ndx(1,cols), 1];
    ijkndx = ijkndx + 1;
    
    % Bottom-left corner
%     winNdx = ndx(rows-1:rows,1:2);
%     i = smallndx(rows/2,1);
%     ijk(ijkndx:ijkndx+4-1,:) = [repmat([i],[4,1]), winNdx(:), reshape(downkern(1:2,2:3)/sum(sum(downkern(1:2,2:3))),[4,1])];
%     ijkndx = ijkndx + 4;
    ijk(ijkndx,:) = [smallndx(floor(rows/2),1), ndx(rows,1), 1];
    ijkndx = ijkndx + 1;
    
    % Bottom-right corner
%     winNdx = ndx(rows-1:rows,cols-1:cols);
%     i = smallndx(rows/2,cols/2);
%     ijk(ijkndx:ijkndx+4-1,:) = [repmat([i],[4,1]), winNdx(:), reshape(downkern(1:2,1:2)/sum(sum(downkern(1:2,1:2))),[4,1])];
%     ijkndx = ijkndx + 4;
    ijk(ijkndx,:) = [smallndx(floor(rows/2),floor(cols/2)), ndx(rows,cols), 1];
    ijkndx = ijkndx + 1;
    
    % Left edge
    for v=3:2:rows-2
%         winNdx = ndx(v-1:v+1,1:2);
%         i = smallndx((v+1)/2,1);
%         ijk(ijkndx:ijkndx+6-1,:) = [repmat([i],[6,1]), winNdx(:), reshape(downkern(:,2:3)/sum(sum(downkern(:,2:3))),[6,1])];
%         ijkndx = ijkndx + 6;
        ijk(ijkndx,:) = [smallndx((v+1)/2,1), ndx(v,1), 1];
        ijkndx = ijkndx + 1;
    end
    
    % Right edge
    for v=3:2:rows-2
%         winNdx = ndx(v-1:v+1,cols-1:cols);
%         i = smallndx((v+1)/2,cols/2);
%         ijk(ijkndx:ijkndx+6-1,:) = [repmat([i],[6,1]), winNdx(:), reshape(downkern(:,1:2)/sum(sum(downkern(:,1:2))),[6,1])];
%         ijkndx = ijkndx + 6;
        ijk(ijkndx,:) = [smallndx((v+1)/2,floor(cols/2)), ndx(v,cols), 1];
        ijkndx = ijkndx + 1;
    end
    
    % Top edge
    for u=3:2:cols-2
%         winNdx = ndx(1:2,u-1:u+1);
%         i = smallndx(1,(u+1)/2);
%         ijk(ijkndx:ijkndx+6-1,:) = [repmat([i],[6,1]), winNdx(:), reshape(downkern(2:3,:)/sum(sum(downkern(2:3,:))),[6,1])];
%         ijkndx = ijkndx + 6;
        ijk(ijkndx,:) = [smallndx(1,(u+1)/2), ndx(1,u), 1];
        ijkndx = ijkndx + 1;
    end
    
    % Bottom edge
    for u=3:2:cols-2
%         winNdx = ndx(rows-1:rows,u-1:u+1);
%         i = smallndx(rows/2,(u+1)/2);
%         ijk(ijkndx:ijkndx+6-1,:) = [repmat([i],[6,1]), winNdx(:), reshape(downkern(1:2,:)/sum(sum(downkern(1:2,:))),[6,1])];
%         ijkndx = ijkndx + 6;
        ijk(ijkndx,:) = [smallndx(floor(rows/2),(u+1)/2), ndx(rows,u), 1];
        ijkndx = ijkndx + 1;
    end

    % Interior pixels
    for u=3:2:cols-1
        for v=3:2:rows-1
            winNdx = ndx(v-1:v+1,u-1:u+1);
            i = smallndx((v+1)/2, (u+1)/2);
            
            downkern = abs(full(A(ndx(v,u),winNdx)));
            %downkern(5) = sum(downkern([1,2,3,4,6,7,8,9]));
            downkern(5) = 0;
            downkern = downkern/sum(downkern(:));
            
            ijk(ijkndx:ijkndx+9-1,:) = [repmat([i],[9,1]), winNdx(:), downkern(:)];
            ijkndx = ijkndx + 9;
        end
    end
    
    down{level} = sparse(ijk(1:ijkndx-1,1), ijk(1:ijkndx-1,2), ijk(1:ijkndx-1,3), (floor(rows/2)*floor(cols/2)), N);
    up{level} = down{level}';
    
    %if( mod(cols,2) == 0 && mod(rows,2) == 0 )
        cols = floor(cols/2);
        rows = floor(rows/2);
    %else
    %    printf('Error');
    %    return;
    %end
    
    A = down{level}*A*up{level};
end

end

