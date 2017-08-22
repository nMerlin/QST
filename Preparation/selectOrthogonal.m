function [O1,O2,O3,oTheta] = selectOrthogonal(X1,X2,X3,theta,piezoSign,varargin)
%SELECTORTHOGONAL Select all data points where X1 and X2 are orthogonal
%
% Input Arguments:
%   (X1,X2,X3) is a 3-Channel dataset, prepared by PREPARE3CHDATA (i.e.
%   reshaped into piezo-segments with offset already removed)
%   theta - reconstructed phase for X3, prepared with COMPUTEPHASE
%   piezoSign - starting sign of the piezo movement
%
% Optional Arguments:
%   selectOrthogonal(...,'Plot','plot'): Visualize the selection process
%   selectOrthogonal(...,'Width',width): Specify the range, in percent
%       of the total data range, to select as orthogonal data points.

%% Constants
SHIFT = 10000;

%% Validate and parse input arguments
p = inputParser;
defaultPlot = 'noplot';
defaultWidth = 0.05; % Width of orthogonal range in % of max-min
addParameter(p,'Plot',defaultPlot,@isstr);
addParameter(p,'Width',defaultWidth,@isnumeric);
parse(p,varargin{:});
c = struct2cell(p.Results);
[plotArg,ORTH_WIDTH] = c{:};

%% Calculate smoothed cross-correlation
ys = smoothCrossCorr(X1,X2);

%% Select triples where X1 and X2 are orthogonal
% Find points, where the cross-correlation of X1 and X2 is in a ORTH_WIDTH
% range around 0. When doing this, we are not distinguishing between a
% phase of pi/2 and a phase of 3*pi/2 between X1 and X2. Therfore, it is
% important to invert the sign of X2 values that are on negative flanks.
globMax = max(ys(:));
globMin = min(ys(:));
halfWidth = ORTH_WIDTH * (globMax - globMin) / 2;
iOrth = find(ys<halfWidth & ys>-halfWidth); % Indices for selection
O1 = X1(iOrth); O2 = X2(iOrth); O3 = X3(iOrth); oTheta = theta(iOrth);
ysShift = circshift(ys,SHIFT,1);
iMinus = find((iOrth>SHIFT) & ((ys(iOrth)-ysShift(iOrth))<0));
ysShift = circshift(ys,-SHIFT,1);
iMinus = [iMinus; find((iOrth<=SHIFT) & ((ys(iOrth)-ysShift(iOrth))>0))];
O2(iMinus) = -O2(iMinus);
O2 = piezoSign*O2;  

%% Visualize the selection process for two piezo segments
if strcmp(plotArg,'plot')
    [nPoints,~] = size(ys);
    x = 1:2*nPoints;
    iOrth = iOrth(iOrth<2*nPoints);
    ys = ys(x);
    plot(x,ys); hold on;
    plot(x(iOrth),ys(iOrth),'.'); hold off;
    legend('Smoothed Cross-Correlation X1.*X2','Selected Data Points');
end

end
