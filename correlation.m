function [ correlationvalue, X ] = correlation( j, data, locs, window )
% CORRELATION calculates <X_i X_i+j>

% eliminate locations where the window would go outside the boundaries
if (locs(1)<ceil(window/2))
    locs = locs(2:end);
elseif ((length(data)-locs(end))<ceil(window/2))
    locs = Locs(1:length(locs)-1);
end

datasize = size(data);
numberOfWindows = length(locs);
X = zeros(numberOfWindows,datasize(2));
correlationvalue = 0;

for n=1:datasize(2)
    % quadrature vector
    points = data(:,n);
    points = points - mean(points);
    for i=1:numberOfWindows
        X(i,n) = x(i, points, locs, window);
    end

    % correlations
    i = 1;
    while (i+j)<=numberOfWindows
        correlationvalue = correlationvalue + X(i,n)*X(i+j,n);
        i = i + 1;
    end
end
correlationvalue = correlationvalue/((numberOfWindows - j)*n);

end

