function [ optimalWindowSize ] = optimalIntegrationWindow( data, locs, maxWindowSize, minWindowSize, stepSize, parallel )
% OPTIMALINTEGRATIONWINDOW computes the correlationvalues <X_i X_i+1> of
% DATA at the locations LOCS, plots them and outputs the minimum
% correlation.

switch nargin
    case 3
        stepSize = 5;
        minWindowSize = 100;
        parallel = 0;
    case 4
        stepSize = 5;
        parallel = 0;
    case 5
        parallel = 0;
end

windowSizes = minWindowSize:stepSize:maxWindowSize;
correlationValues = zeros(length(windowSizes),1);

numberOfWindowSizes = length(windowSizes);

if parallel==1
    parfor i=1:numberOfWindowSizes
        correlationValues(i) = correlation(1, data, locs, windowSizes(i));
        disp(['Windowsize:',num2str(i),'/',num2str(numberOfWindowSizes)]);
    end
else
    for i=1:numberOfWindowSizes
        correlationValues(i) = correlation(1, data, locs, windowSizes(i));
        disp(['Windowsize:',num2str(i),'/',num2str(numberOfWindowSizes)]);
    end
end

plot(windowSizes,correlationValues);
[~, index] = min(abs(correlationValues));
optimalWindowSize = windowSizes(index);
txtOptimalWindowSize = strcat('\leftarrow', num2str(optimalWindowSize));
text(optimalWindowSize,correlationValues(index),txtOptimalWindowSize);

end
