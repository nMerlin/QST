function [ optimalWindowSize ] = optimalIntegrationWindow(plotFlag, data, locs, maxWindowSize, minWindowSize, stepSize, parallel )
% OPTIMALINTEGRATIONWINDOW computes the correlationvalues <X_i X_i+1> of
% DATA at the locations LOCS for different window sizes, plots them and
% outputs the minimum correlation.
%
%   Special Inputs: If PARALLEL=1, the function uses PARFOR instead of FOR.
%   STEPSIZE gives the stepwidth used for the plot and calculation. If
%   PLOTFLAG=1, the function also plots the correlationvalues and writes
%   them to the file optimalIntegrationWindow.png
%
%   See also: PARFOR

switch nargin
    case 4
        stepSize = 5;
        minWindowSize = 200;
        parallel = 0;
    case 5
        stepSize = 5;
        parallel = 0;
    case 6
        parallel = 0;
end

windowSizes = minWindowSize:stepSize:maxWindowSize;
correlationValues = zeros(length(windowSizes),1);

numberOfWindowSizes = length(windowSizes);

if parallel==1
    parfor i=1:numberOfWindowSizes
        correlationValues(i) = correlation(1, data, locs, windowSizes(i));
        disp(['Windowsize:',num2str(i/numberOfWindowSizes),'%']);
    end
else
    for i=1:numberOfWindowSizes
        correlationValues(i) = correlation(1, data, locs, windowSizes(i));
        disp(['Windowsize:',num2str(i/numberOfWindowSizes*100),'%']);
    end
end

[~, index] = min(abs(correlationValues));
optimalWindowSize = windowSizes(index);

if plotFlag==1
    plot(windowSizes,correlationValues);
    txtOptimalWindowSize = strcat('\leftarrow', num2str(optimalWindowSize));
    text(optimalWindowSize,correlationValues(index),txtOptimalWindowSize);
    print('optimalIntegrationWindow.png','-dpng');
end

end
