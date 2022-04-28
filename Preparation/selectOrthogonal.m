function [O1,O2,O3,oTheta,iOrth,ys] = selectOrthogonal(X1,X2,X3,theta,piezoSign,varargin)
%SELECTORTHOGONAL Select all data points where X1 and X2 are orthogonal
%
% Input Arguments:
%   (X1,X2,X3) is a 3-Channel dataset, prepared by PREPARE3CHDATA (i.e.
%   reshaped into piezo-segments with offset already removed)
%   X1 is the channel that is used for orthogonal selection as well as for 
%   the relative phase computation between X1 and X3. 
%   X2 is only used for orthogonal selection. X3 is the target channel.
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
defaultWidth = 0.05; % Width of orthogonal range as ratio of max-min
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
% phase of pi/2 and a phase of 3*pi/2 between X1 and X2. Therefore, it is
% important to invert the sign of some X2 values: Either of those values on
% negative flanks or positive ones. In this case, we invert the values on
% the negative flanks.
globMax = max(ys(:));
globMin = min(ys(:));
halfWidth = ORTH_WIDTH * (globMax - globMin) / 2;
iOrth = find(ys<halfWidth & ys>-halfWidth); % Indices for selection

%%besser machen
% [nbunch,nSeg]=size(ys);
% X1 = reshape(X1,[nbunch,nSeg]);X2 = reshape(X2,[nbunch,nSeg]);
% [O1,O2,iOrth]=deal(NaN([nbunch,nSeg]));
% for iSeg = 1:nSeg
%     locMax = max(ys(:,iSeg));
%     locMin = min(ys(:,iSeg));
%     halfWidth = ORTH_WIDTH * (locMax - locMin) / 2;
%     iOrthSeg = find(ys(:,iSeg)<halfWidth & ys(:,iSeg)>-halfWidth); % Indices for selection
%     O1(1:length(iOrthSeg),iSeg)=X1(iOrthSeg,iSeg);
%     O2(1:length(iOrthSeg),iSeg)=X2(iOrthSeg,iSeg);
%     iOrth(1:length(iOrthSeg),iSeg)=iOrthSeg;%geht nicht
% end
% O1=O1(:);O2=O2(:);iOrth=iOrth(:);
% O1=O1(~isnan(O1));O2=O2(~isnan(O2));iOrth=iOrth(~isnan(iOrth));
%%besser machen

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
    x = 1:238*nPoints;
    iOrthPlot = iOrth(iOrth<238*nPoints);
    ysPlot = ys(x);
    plot(x,ysPlot); hold on;
    plot(x(iOrthPlot),ysPlot(iOrthPlot),'.'); hold off;
    legend('Smoothed Cross-Correlation X1.*X2','Selected Data Points');
end

end
