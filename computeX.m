function [ x_i ] = computeX( i, points, locs, window )
% COMPUTEX adds up the data points in POINTS lying in the window with size
% WINDOW symmetrically around the I-th location in LOCS and divides the
% results by the number of points in the integration window.
%
% Inputs: POINTS is an array containing doubles, LOCS is an array
% containing integers and WINDOW is an integer value.

rows = length(points);

start = locs(i)-ceil(window/2);
stop = locs(i)+ceil(window/2);

assert(start>0,'Start of window is not positive!');
assert(stop<(rows+1),'Stop of window is out of range!');
x_i = sum(points(start:stop))/(stop-start+1);

end

