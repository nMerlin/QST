function plotSeriesPostselections(listOfParams,varargin)
%PLOTSERIESPOSTSELECTIONS

%% Validate and parse input arguments
p = inputParser;
defaultFilename = '';
addParameter(p,'Filename',defaultFilename,@isstr);
defaultType = 'Amplitude';
addParameter(p,'Type',defaultType,@isstr);
parse(p,varargin{:});
c = struct2cell(p.Results);
[filename,typestr] = c{:};

%% Gather data
[X,Y,Z] = deal([]);
for iParams = 1:length(listOfParams)
    selParams = listOfParams(iParams);
    A = seriesRead3ChTable(selParams);
    H = height(A);
    Radii = ones(H,1) * selParams.Position(1);
    X(iParams,:) = A.Delay;
    [X(iParams,:),I] = sort(X(iParams,:)); % Sort for Delays
    Y(iParams,:) = Radii;
    Z(iParams,:) = A.discAmpl;
    Z(iParams,:) = Z(iParams,I);
end
[~,I] = sort(X(:,1)); % Sort for Radii
X = X(I,:);
Y = Y(I,:);
Z = Z(I,:);

%% Create figure
fig = figure;
formatFigA5(fig);
switch typestr
    case 'Amplitude'
        waterfall(X,Y,Z);
        view(-20,20);
end
set(fig,'Color','w');

%% Write figure to file and close it
if ~isempty(filename)
    savefig(fig,filename);
    close all;
end

end

