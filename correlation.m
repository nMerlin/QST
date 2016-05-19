function [ correlationvalue, X ] = correlation( j, data, locs, window )
% CORRELATION calculates <X_i X_i+j>_i for the data arrays given in the
% columns of DATA.
%
% Outpus: CORRELATIONVALUE is the resulting correlation <X_i X_i+j>_i,
% averaged over all arrays and all locations i in LOCS. X is a matrix
% containing the extracted normalized sums for all suitable windows in
% DATA.
%
% See also: POINTWISEVARIANCE

% Eliminate locations whose corresponding window would be outside the range
% of DATA.
if (locs(1)<ceil(window/2))
    locs = locs(2:end);
elseif ((length(data)-locs(end))<ceil(window/2))
    locs = locs(1:length(locs)-1);
end

[~,columns] = size(data);
numberOfWindows = length(locs);
X = zeros(numberOfWindows,columns);
correlationvalue = 0;

for column=1:columns
    % Compute quadrature vector
    points = data(:,column);
    points = points - mean(points);
    for i=1:numberOfWindows
        X(i,column) = computeX(i, points, locs, window);
    end

    % Compute CORRELATIONVALUE
    i = 1;
    while (i+j)<=numberOfWindows
        correlationvalue = correlationvalue + X(i,column)*X(i+j,column);
        i = i + 1;
    end
end

% Normalization
correlationvalue = correlationvalue/((numberOfWindows - j)*columns);

end

