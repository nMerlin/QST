function [X,theta,iSelect] = selectRegionOfTotalPhase(O1,O2,O3,theta,thetaMira,varargin)
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
    case 'phase'
        mTheta = position(1);
        dTheta = position(2);
        husimiPhase = pi + atan2(O2,O1);
        thetaTotal = mod(husimiPhase + theta - thetaMira,2*pi);
        iSelect = find(thetaTotal >mTheta-dTheta/2 & ...
            thetaTotal <mTheta+dTheta/2 & (O1.^2+O2.^2)>0.5);
     case 'phaseAndAmplitude'
        mTheta = position(1);
        dTheta = position(2);
        r = position(3);
        w = position(4);
        husimiPhase = pi + atan2(O2,O1);
        thetaTotal = mod(husimiPhase + theta - thetaMira,2*pi);
        iSelect = find(thetaTotal>mTheta-dTheta/2 & ...
            thetaTotal<mTheta+dTheta/2 & (O1.^2+O2.^2)>0.5 & ...
            sqrt(O1.^2+O2.^2)>r-w/2 & sqrt(O1.^2+O2.^2)<r+w/2);
end
X = O3(iSelect);

%% Phase calculation

theta = theta(iSelect);

%% Show selection
if strcmp(plotopt,'show')
    assessTheta(theta,X,'Husimi',{O1,O2,iSelect},'VarBins', ...
        200,'PhaseBins',200,'Output',output,'Filename',filename);
end

end
