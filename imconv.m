function [ out ] = imconv( I, kernel )
%IMCONV Summary of this function goes here
%   Detailed explanation goes here

out = cat(3, conv2(I(:,:,1),kernel,'same'), conv2(I(:,:,2),kernel,'same'), conv2(I(:,:,3),kernel,'same') );

end

