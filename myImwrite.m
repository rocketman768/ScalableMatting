function [ ] = myImwrite( I, filename, range, map )
%MYIMWRITE Write 2D matrix to image file.

Imin = min(min(I));
Imax = max(max(I));

% Scale to [0,1]
if( ~exist('range','var') || isempty(range) )
    J = double(I-Imin)/(Imax-Imin);
else
    J = double(I-range(1))/(range(2)-range(1));
end

if( ~exist('map','var') )
    imwrite(J,filename);
else
    imwrite(1+J*(size(map,1)-1), map, filename);
end

end

