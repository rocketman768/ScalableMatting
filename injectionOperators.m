function [ up, down ] = injectionOperators( I, numLevels )
%INJECTIONOPERATORS Summary of this function goes here
%   Detailed explanation goes here

rows = size(I,1);
cols = size(I,2);

for level=numLevels:-1:1
    ndx = reshape(1:rows*cols,[rows,cols]);
    rowNdx = (1:rows)';
    colNdx = (1:cols);
    selected = logical(bitand(rowNdx,1) * bitand(colNdx,1));
    numSelected = sum(selected(:));
    
    down{level} = sparse(1:numSelected, ndx(selected(:)), ones(numSelected,1), numSelected, rows*cols);
    up{level} = down{level}';
    
    rows = sum(bitand(rowNdx,1));
    cols = sum(bitand(colNdx,1));
end

end

