function [ frame ] = myIm2Frame( im, map )
%MYIM2FRAME Convert grayscale 0 < im < 1 to a movie frame.

N = size(map,1);

im(im < 0) = 0;
im(im > 1) = 1;

frame = im2frame( im*(N-1)+1, map );

end

