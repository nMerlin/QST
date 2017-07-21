function [X,theta] = selectRegion(O1,O2,O3,theta,varargin)
%SELECTREGION Return all values of O3 in the specified (O1,O2)-region
%
%   Rectangle:
%   selectRegion(O1,O2,O3,theta,'rectangle',[x y w h]) - lower left corner
%   at (x,y) and (width,height) = (w,h)

%% Validate and parse input arguments
p = inputParser;
defaultType = 'rectangle';
defaultPosition = [2 2 0.5 0.5];
defaultPlotOpt = 'hide';
addParameter(p,'Type',defaultType,@isstr);
addParameter(p,'Position',defaultPosition,@isvector);
addParameter(p,'Plot',defaultPlotOpt,@isstr);
parse(p,varargin{:});
c = struct2cell(p.Results);
[plotopt,position,type] = c{:};

%% Selection process
switch type
    case 'rectangle'
        [x,y,w,h] = position(1:4);
        iSelect = find(O1>x & O1<x+w & O2>y & ...
            O2<y+h);
    case 'halfcircle'
        r = position(1); % inner radius
        w = position(2); % outer radius = r+w
        iSelect = find(sqrt(O1.^2+O2.^2)>r & sqrt(O1.^2+O2.^2)<r+w & O2>0);
end
X = O3(iSelect);

%% Phase calculation
quot = O2./O1;
iTan = find(O1>0);
theta(iTan) = theta(iTan) - (atan(quot(iTan)));
iTan = find(O1<=0);
theta(iTan) = theta(iTan) - (atan(quot(iTan)) + pi);
theta = theta(iSelect);

%% Show selection
if strcmp(plotopt,'show')
    [H, binsO1, binsO2] = histogram2D(O1,O2);
    imagesc(binsO1,binsO2,H); axis on; colormap hot; hold on;
    plot(O1(iSelect),O2(iSelect),'.'); hold off;
end

end
