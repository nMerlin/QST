function [] = plotFromPng(filename)
A = imread(filename);
image(A);
end