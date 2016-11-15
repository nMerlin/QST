function [ correlationvalue, X ] = correlation( j, data, locs )
% CORRELATION calculates <X_i X_i+j>_i for the data arrays given in the
% columns of DATA.
%
% Outputs: CORRELATIONVALUE is the resulting correlation <X_i X_i+j>_i,
% averaged over all arrays and all locations i in LOCS. X is a matrix
% containing the extracted normalized sums for all suitable windows in
% DATA.
%
% See also: POINTWISEVARIANCE

WINDOWSIZE = 2/3; % Integrate over WINDOWSIZE * mean(diff(locs))

% Calculating the size of the integration window

% Eliminate locations whose corresponding window would be outside the range
% of DATA.
window = round(WINDOWSIZE * mean(diff(locs)));
if (locs(1)<ceil(window/2))
    locs = locs(2:end);
end
if ((length(data)-locs(end))<ceil(window/2))
    locs = locs(1:length(locs)-1);
end

[~,columns] = size(data);
numberOfWindows = length(locs);
X = zeros(numberOfWindows,columns);
correlationvalues = zeros(2,2);

for column=1:columns
    % Compute quadrature vector
    points = data(:,column);
    %points = points - mean(points);
    for i=1:numberOfWindows
        X(i,column) = computeX(i, points, locs, window);
    end

    % Compute CORRELATIONVALUE
    i = 1;
    A = zeros(1,numberOfWindows-j);
    B = A;
    while (i+j)<=numberOfWindows
        A(i)=X(i,column);
        B(i)=X(i+j,column);
        i = i + 1;
    end
    correlationvalues = correlationvalues + corrcoef(A,B);
end

% Normalization
correlationvalues = correlationvalues./columns;
correlationvalue = correlationvalues(1,2);

end

