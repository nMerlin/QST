function [X,theta] = selectRegion(O1,O2,O3,theta,varargin)
%SELECTREGION Return all values of O3 in the specified (O1,O2)-region
%
%   Rectangle:
%   selectRegion(O1,O2,O3,theta,'rectangle',[x y w h]) - lower left corner
%   at (x,y) and (width,height) = (w,h)

%% Handle optional input arguments and default values
nVarargin = length(varargin);
optArgs = {'rectangle' 2 2 0.5 0.5};
optArgs(1:nVarargin) = varargin;
[type,x,y,w,h] = optArgs{:}; % Smoothing parameter

%% Selection process
switch type
    case 'rectangle'
        iSelect = find(O1>x & O1<x+w & O2>y & ...
            O2<y+h);
end
X = O3(iSelect);

%% Phase calculation
quot = O2./O1;
iTan = find(O1>0);
theta(iTan) = theta(iTan) + atan(quot(iTan));
iTan = find(O1<=0);
theta(iTan) = theta(iTan) + atan(quot(iTan)) + pi;
theta = theta(iSelect);

end
