function g2 = g2( X, nResolution )
%G2 Creates a plot showing the g2(0) behavior over time
%
% Input Parameters:
%   X - Quadratures from a continuous quantum state measurement (equal time
%       spacing between all points)
%   NRESOLUTION - Number of Quadratures used to create a single data point
%           (time resolution)
%   FILENAME - Output filename for the resulting plot

%%% Reshaping X according to NRESOLUTION
X = X(:);
nSegments = floor(length(X)/nResolution);
X = X(1:nSegments*nResolution);
X = reshape(X,[nResolution nSegments]); % Consecutive values in columns

%%% Piecewise photon number <a^+ a>
X = X - mean(mean(X));
ada = mean(X.^2)-0.5;

%%% Piecewise <a^+ a^+ a a>, in the following adadaa, and g2(t,t)
adadaa = 2/3*mean(X.^4)-2*ada-0.5;
g2 = adadaa./ada.^2;
g2 = g2';

end

