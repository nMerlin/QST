function [ fitParams, fitFunction,fval, exitFlag] = fitSinus( x, y, varargin)
%fitSinusodial Fits a sine function to a (x,y) dataset
%
%   (x,y) should already be sinusodial for the fit to work properly.
%   varargin = 'rmLin': additionally account for a linear trend
%              'show': show plot of raw data and fitted sinusoidal function
%
%   FITPARAMS outputs the fit paramters b of b1*sin(x)+b4
%   FITFUNCTION is the function handle of the above fit function
%   EXITFLAG is the exit status of FMINSEARCH used to find optimal
%   parameters
%
% Adapted from Star Strider:
% https://www.mathworks.com/matlabcentral/answers/ ...
% 121579-curve-fitting-to-a-sinusoidal-function#answer_128539

    % Optional Parameters
    rmLin = 0;
    show = 0;
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
    yOffset = yMin + yRange/2; % Estimate offset

    % Function to fit
    fitFunction = @(b,x)  b(1).*(sin(x)) + b(2);
    % Least-Squares cost function
    squaresFunction = @(b) sum((fitFunction(b,x) - y).^2);
    % Minimise Least-Squares (with lower boundaries)
    [fitParams, fval, exitFlag] = fminsearchbnd(squaresFunction, ...
        [yRange/2; yOffset], [0; -inf]);

    % Correct for linear trend
    if rmLin == 1
        fitFunction = @(b,x)  b(1).*(sin(x )) ...
            + b(2) + b(3) .* x;
        squaresFunction = @(b) sum((fitFunction(b,x) - y).^2);
        [fitParams, fval, exitFlag] = fminsearchbnd(squaresFunction, ...
        [fitParams; 0], [0;  -inf; -inf]);
    end
    
    %%% Plot raw data and fitted function
    if show == 1
        figure(1)
        plot(x,y,'b',x,fitFunction(fitParams,x),'r')
        grid
    end
end

