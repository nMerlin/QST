function [X,theta,iSelect] = selectRegion(O1,O2,O3,theta,varargin)
%SELECTREGION Return all values of O3 in the specified (O1,O2)-region
%
%   Rectangle:
%   selectRegion(O1,O2,O3,theta,'rectangle',[x y w h]) - lower left corner
%   at (x,y) and (width,height) = (w,h)
%
% Options:
%   Output - 'print' saves plotted figure to file

%% Validate and parse input arguments
p = inputParser;
defaultType = 'dot';
defaultPosition = [2 2 0.25];
defaultPlotOpt = 'hide';
defaultOutput = 'figure';
defaultFilename = '';
addParameter(p,'Type',defaultType,@isstr);
addParameter(p,'Position',defaultPosition,@isvector);
addParameter(p,'Plot',defaultPlotOpt,@isstr);
addParameter(p,'Output',defaultOutput,@isstr);
addParameter(p,'Filename',defaultFilename,@isstr);
parse(p,varargin{:});
c = struct2cell(p.Results);
[filename,output,plotopt,position,type] = c{:};

%% Selection process
switch type
    case 'rectangle'
        x = position(1);
        y = position(2);
        w = position(3);
        h = position(4);
        iSelect = find(O1>x & O1<x+w & O2>y & ...
            O2<y+h);
    case 'halfcircle'
        r = position(1); % inner radius
        w = position(2); % outer radius = r+w
        iSelect = find(sqrt(O1.^2+O2.^2)>r & sqrt(O1.^2+O2.^2)<r+w & O2>0);
    case 'fullcircle'
        r = position(1); % central radius
        w = position(2); % circle thickness
        iSelect = find(sqrt(O1.^2+O2.^2)>r-w/2 & sqrt(O1.^2+O2.^2)<r+w/2);
    case 'dot'
        x = position(1);
        y = position(2);
        r = position(3);
        iSelect = find(((O1-x).^2+(O2-y).^2) <= r^2);
    case 'phase'
        mTheta = position(1);
        dTheta = position(2);
        husimiPhase = atan2(O2,O1);
        iSelect = find(husimiPhase>mTheta-dTheta/2 & ...
            husimiPhase<mTheta+dTheta/2 & (O1.^2+O2.^2)>0.5);
    case 'phaseAndAmplitude'
        mTheta = position(1);
        dTheta = position(2);
        r = position(3);
        w = position(4);
        husimiPhase = atan2(O2,O1);
        iSelect = find(husimiPhase>mTheta-dTheta/2 & ...
            husimiPhase<mTheta+dTheta/2 & (O1.^2+O2.^2)>0.5 & ...
            sqrt(O1.^2+O2.^2)>r-w/2 & sqrt(O1.^2+O2.^2)<r+w/2);
    case 'Qline'
        Q = position(1); % line center
        w = position(2); % line width
        iSelect = find(O1>Q-w/2 & O1<Q+w/2);
    case 'Pline'
        P = position(1); % line center
        w = position(2); % line width
        iSelect = find(O2>P-w/2 & O2<P+w/2);
end
X = O3(iSelect);

%% Phase calculation
theta = mod(theta + atan2(O2,O1),2*pi);
theta = theta(iSelect);

%% Show selection
if strcmp(plotopt,'show')
    assessTheta(theta,X,'Husimi',{O1,O2,iSelect},'VarBins', ...
        200,'PhaseBins',200,'Output',output,'Filename',filename);
end

end
