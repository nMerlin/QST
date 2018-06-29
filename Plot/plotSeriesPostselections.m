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
[X,Y,discAmpl,discMeanVar,discN] = deal([]);
for iParams = 1:length(listOfParams)
    selParams = listOfParams(iParams);
    A = seriesRead3ChTable(selParams);
    H = height(A);
    Radii = ones(H,1) * selParams.Position(1);
    X(iParams,:) = A.Delay;
    [X(iParams,:),I] = sort(X(iParams,:)); % Sort for Delays
    Y(iParams,:) = Radii;
    discAmpl(iParams,:) = A.discAmpl;
    discAmpl(iParams,:) = discAmpl(iParams,I);
    discMeanVar(iParams,:) = A.discMeanVar;
    discMeanVar(iParams,:) = discMeanVar(iParams,I);
    discN(iParams,:) = A.discN;
    discN(iParams,:) = discN(iParams,I);
end
[~,I] = sort(X(:,1)); % Sort for Radii
X = X(I,:);
Y = Y(I,:);
discAmpl = discAmpl(I,:);
discMeanVar = discMeanVar(I,:);
discN = discN(I,:);

%% Create figure
fig = figure;
formatFigA5(fig);
switch typestr
    case 'Amplitude'
        waterfall(X,Y,discAmpl);
        view(-20,20);
        xlabel('Delay (fs)');
        zlabel('Coherent Amplitude');
        title('Coherent Amplitude vs. Radius of Postselected Fullcircle');
    case 'MeanVar'
        surf(X,Y,discMeanVar);
        view(-50,20);
        xlabel('Delay (fs)');
        zlabel('Average Variance');
        title('Variance vs. Radius of Postselected Fullcircle');
    case 'DiscN'
        waterfall(X,Y,discN);
        view(-20,20);
        xlabel('Delay (fs)');
        zlabel('Photon Number');
        title('Photon Number vs. Radius of Postselected Fullcircle');
end
set(fig,'Color','w');

%% Write figure to file and close it
if ~isempty(filename)
    savefig(fig,filename);
    close all;
end

end

