function [X,theta] = selectRegion(O1,O2,O3,theta,varargin)
%SELECTREGION Return all values of O3 in the specified (O1,O2)-region
%
%   Rectangle:
%   selectRegion(O1,O2,O3,theta,'rectangle',[x y w h]) - lower left corner
%   at (x,y) and (width,height) = (w,h)

%% Handle optional input arguments and default values
nVarargin = length(varargin);
optArgs = {'rectangle' [2 2 0.5 0.5]};
optArgs(1:nVarargin) = varargin;
[type pos] = optArgs{:}; % Smoothing parameter

%% Selection process
switch type
    case 'rectangle'
        iSelect = find(O1>pos(1) & O1<pos(1)+pos(3) & O2>pos(2) & ...
            O2<pos(2)+pos(4));
end
X = O3(iSelect);
theta = theta(iSelect);

end
