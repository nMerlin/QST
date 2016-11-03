function [ fitParams, fitFunction, exitFlag] = fitSinusodial( x, y )
%fitSinusodial Fits a sine function to a (x,y) dataset
%   (x,y) should already be sinusodial for the fit to work properly.
%   FITPARAMS outputs the fit paramters b of b1*sin(2*pi/b2*x+2*pi/b3)+b4
%   FITFUNCTION is the function handle of the above fit function
%   EXITFLAG is the exit status of FMINSEARCH used to find an optimal
%   parameters
%
% Adapted from Star Strider:
% https://www.mathworks.com/matlabcentral/answers/ ...
% 121579-curve-fitting-to-a-sinusoidal-function#answer_128539

    y = y(~isnan(y)); % Remove NaN-entries in y
    x = x(~isnan(x));
    if ~isrow(y)
        y = y';
    end
    if ~isrow(x)
        x = x';
    end
    
    yMax = max(y);
    yMin = min(y);
    yRange = (yMax-yMin);
    yZero = y-yMax+(yRange/2);
    zeroIndicator = yZero.* circshift(yZero, [0 1]); ...
        % Find zero-crossings (ignore first entry in zeroIndicator)
    zerosX = x([1 zeroIndicator(2:end)] <= 0);
    period = 2*mean(diff(zerosX));
    yOffset = yMin + yRange/2; % Estimate offset
    yPhase = mod(-period/zerosX(1),period); % Estimate phase

    % Function to fit
    fitFunction = @(b,x)  b(1).*(sin(2*pi*x./b(2) + 2*pi/b(3))) + b(4);
    % Least-Squares cost function
    squaresFunction = @(b) sum((fitFunction(b,x) - y).^2);
    % Minimise Least-Squares
    [fitParams, ~, exitFlag] = fminsearch(squaresFunction, ...
        [yRange/2; period;  yPhase;  yOffset]);

    %%% Plot
    figure(1)
    plot(x,y,'b',x,fitFunction(fitParams,x),'r')
    grid
end

