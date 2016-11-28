function [ fitParams, fitFunction, exitFlag] = fitSinusoidal( x, y, varargin)
%fitSinusodial Fits a sine function to a (x,y) dataset
%
%   (x,y) should already be sinusodial for the fit to work properly.
%   varargin = 'rmLin': additionally account for a linear trend
%
%   FITPARAMS outputs the fit paramters b of b1*sin(2*pi/b2*x+2*pi/b3)+b4
%   FITFUNCTION is the function handle of the above fit function
%   EXITFLAG is the exit status of FMINSEARCH used to find optimal
%   parameters
%
% Adapted from Star Strider:
% https://www.mathworks.com/matlabcentral/answers/ ...
% 121579-curve-fitting-to-a-sinusoidal-function#answer_128539

    % Optional Parameters
    rmLin = 0;
    if nargin > 2
        for i = 3:nargin
            eval([varargin{i-2} '=1;']);
        end
    end

    % Correcting input data format
    y = y(~isnan(y)); % Remove NaN-entries in y
    x = x(~isnan(x));
    if ~isrow(y)
        y = y';
    end
    if ~isrow(x)
        x = x';
    end
    
    % Estimate fit paramters
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

    % Correct for linear trend, if applicable
    if rmLin == 1
        fitFunction = @(b,x)  b(1).*(sin(2*pi*x./b(2) + 2*pi/b(3))) ...
            + b(4) + b(5) .* x;
        squaresFunction = @(b) sum((fitFunction(b,x) - y).^2);
        [fitParams, ~, exitFlag] = fminsearch(squaresFunction, ...
        [fitParams; 0]);
    end
    
    %%% Plot
    figure(1)
    plot(x,y,'b',x,fitFunction(fitParams,x),'r')
    grid
end

