function [ x0 ] = initialGuess( scribs )
%INITIALGUESS Make a good initial guess for the solution given the scribs
%   scribs - scribble image

% x0(i) is the value of the closest constraint.
scribs = scribs(:,:,1);
[D,idx] = bwdist(scribs > 0.9 | scribs < 0.1, 'euclidean');
x0 = scribs(idx);
x0 = x0(:);

end

