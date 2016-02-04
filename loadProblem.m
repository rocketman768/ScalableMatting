% Stuff for user to set++++++++++++++++++++++++++++++++++++++++++++++++++++
I = double(imread('highres/4Mpx/27.ppm'))./255;
alphagt = double(imread('highres/4Mpx/27_gt.pgm'))./255;
scribs = double(imread('highres/4Mpx/27_scribs.pgm'))./255;

eps = 1e-3;
gamma = 1;
%--------------------------------------------------------------------------

alphagt = alphagt(:);

[c,v] = scribData(scribs');
c = c(:);
v = v(:);

b = gamma*v;

N = size(I,1)*size(I,2);

[up,down] = naiveTransferOperators(I,8,true);
