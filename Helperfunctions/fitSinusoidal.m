function [fitParams,fitFunction,exitFlag] = fitSinusoidal(x,y,varargin)
%FITSINUSOIDAL Fits a sine function to a (x,y) dataset. (x,y) should
%already be sinusoidal for the fit to work properly. The dataset will be
%sorted according to the x values before fitting.
%
% Output Parameters:
%   FITPARAMS outputs the fit paramters b of b1*sin(2*pi/b2*x+2*pi/b3)+b4
%   FITFUNCTION is the function handle of the above fit function
%   EXITFLAG is the exit status of FMINSEARCH used to find optimal
%   parameters
%
% Optional Input Parameters:
%   'RemoveLinearOffset': Default is false. If true, a linear trend is
%       removed from the sine function.
%   'Plot': Default is 'hide'. If 'show', a graphical output is generated.
%   'Periods': Number of periods, default is 1. When -1, the algorithm
%       tries to find the number of periods according to the zero crossings
%       it finds (experimental).
%   'Phase': Phase of first datapoint (in rad), default is 0.
%
% Adapted from Star Strider:
% https://www.mathworks.com/matlabcentral/answers/ ...
% 121579-curve-fitting-to-a-sinusoidal-function#answer_128539

%% Validate and parse input arguments
p = inputParser;
defaultPeriods = 1;
addParameter(p,'Periods',defaultPeriods,@isnumeric);
defaultPhase = 0;
addParameter(p,'Phase',defaultPhase,@isnumeric);
defaultPlotOpt = 'hide';
addParameter(p,'Plot',defaultPlotOpt,@isstr);
defaultRemoveLinearOffset = false;
addParameter(p,'RemoveLinearOffset',defaultRemoveLinearOffset,@islogical);
parse(p,varargin{:});
c = struct2cell(p.Results);
[periods,phase,plotOpt,rmFlag] = c{:};

%% Correcting input data format
% Remove NaN-entries
y = y(~isnan(y));
x = x(~isnan(x));
% Ensure row-vectors
if ~isrow(y)
    y = y';
end
if ~isrow(x)
    x = x';
end
% Sort according to x values
[x,sortedI] = sort(x);
y = y(sortedI);

%% Estimate fit paramters
% Estimate offset
yMax = max(y);
yMin = min(y);
yRange = (yMax-yMin);
yZero = y-yMax+(yRange/2);
yOffset = yMin + yRange/2;
if periods == -1
    % Estimate periods by finding zero crossings
    zeroIndicator = yZero.* circshift(yZero, [0 1]);
    zerosX = x([1 zeroIndicator(2:end)] <= 0);
    period = 2*mean(diff(zerosX));
else
    period = (x(end)-x(1))/periods;
end

% Function to fit
fitFunction = @(b,x)  b(1).*(sin(2*pi*x./b(2) + b(3))) + b(4);
% Least-Squares cost function
squaresFunction = @(b) sum((fitFunction(b,x) - y).^2);
% Minimize Least-Squares (with lower boundaries)
[fitParams,~,exitFlag] = fminsearchbnd(squaresFunction, ...
    [yRange/2;period;phase;yOffset], [0;0;-inf;-inf]);

% Correct for linear trend
if rmFlag
    fitFunction = @(b,x)  b(1).*(sin(2*pi*x./b(2) + 2*pi/b(3))) ...
        + b(4) + b(5) .* x;
    squaresFunction = @(b) sum((fitFunction(b,x) - y).^2);
    [fitParams,~,exitFlag] = fminsearchbnd(squaresFunction, ...
    [fitParams; 0], [0; 0; -inf; -inf; -inf]);
end

%% Plot raw data and fitted function
if strcmp(plotOpt,'show')
    figure(1)
    plot(x,y,'.',x,fitFunction(fitParams,x),'r')
    grid
end
    
end % fitSinusoidal

