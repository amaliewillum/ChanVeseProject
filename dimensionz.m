function [P] = dimensionz(im)
P = imread(im);
if ndims(P) ~= 2
    P = rgb2gray(P);
end
P = double(P);
