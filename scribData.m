function [ constrained, vals ] = scribData( scribImage )

if( ischar(scribImage) )
    scribs = double(imread(scribImage)) ./ 255;
else
    scribs = scribImage;
end
scribs = scribs(:,:,1);
% scribs(scribs<0.1) = 0;
% scribs(scribs>0.9) = 1;
% 
% constrained = (scribs < 0.1);
% constrained = constrained + (scribs > 0.9);

constrained = (scribs == 0);
constrained = constrained + (scribs == 1);

%vals = constrained .* scribs;
vals = constrained .* scribs;
