function [O1,O2,O3] = selectOrthogonal(X1,X2,X3,varargin)
%SELECTORTHOGONAL Select all data points where X1 and X2 are orthogonal
%
% Input Arguments:
%   (X1,X2,X3) is a 3-Channel dataset, prepared by PREPARE3CHDATA (i.e.
%   reshaped into piezo-segments with offset already removed)
%
% Optional Arguments:
%   selectOrthogonal(X1,X2,X3,ORTH_WIDTH): Specify the range, in percent of
%       the total data range, to select as orthogonal data points.

%% Handle optional input arguments and default values
nVarargin = length(varargin);
optArgs = {0.05};
optArgs(1:nVarargin) = varargin;
[ORTH_WIDTH] = optArgs{:}; %Width of orthogonal range in % of max-min

%% Calculate smoothed cross-correlation
ys = smoothCrossCorr(X1,X2);

%% Select triples where X1 and X2 are orthogonal
globMax = max(ys(:));
globMin = min(ys(:));
halfWidth = ORTH_WIDTH * (globMax - globMin) / 2;
iOrth = find(ys<halfWidth & ys>-halfWidth); % Indices for selection
O1 = X1(iOrth); O2 = X2(iOrth); O3 = X3(iOrth);

end

