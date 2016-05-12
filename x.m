function [ x_i ] = x( i, points, locs, window )
%X derives the i-th quadrature value for the given detuning and window-size

start = locs(i)-ceil(window/2);
stop = locs(i)+ceil(window/2);
x_i = sum(points(start:stop));

end

