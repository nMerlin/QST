function makeSelectionPlots(type,varargin)
%MAKESELECTIONPLOTS

%% Validate and parse input arguments
p = inputParser;
defaultListOfParams = '';
addParameter(p,'ListOfParams',defaultListOfParams,@isstr);
parse(p,varargin{:});
c = struct2cell(p.Results);
[listOfParams] = c{:};

% Constants
pdfpath = 'figures-pdf/';
figurepath = 'figures-fig/';

if nargin == 0
    type = 'all';
end
if isempty(listOfParams)
    listMeanVarX = dir('*-series3Ch-*');
    listOfParams = cellfun(@selStrToParams,{listMeanVarX.name});
end

[radiusDiscAmpl,cleanRadiusDiscAmpl,cleanDelayPlots, ...
    delayPlots,pdfs,cleanpdfs] = deal(false);
switch type
    case 'all'
        delayPlots = true;
        radiusDiscAmpl = true;
    case 'radiusPlots'
        radiusDiscAmpl = true;
    case 'cleanRadiusPlots'
        cleanRadiusDiscAmpl = true;
    case 'delayPlots'
        delayPlots = true;
    case 'cleanDelayPlots'
        cleanDelayPlots = true;
    case 'pdfs'
        pdfs = true;
    case 'cleanpdfs'
        cleanpdfs = true;
end

%% Resolve dependencies

%% Make stuff
datestring = datestr(date,'yyyy-mm-dd');
if delayPlots
    for iParams = 1:length(listOfParams)
        makeDelayPlots('plots','SelectionParameters',listOfParams(iParams));
    end
end
if radiusDiscAmpl
    filenameFig = [figurepath,datestring,'-RadiusDiscAmpl.fig'];
    plotSeriesPostselections(listOfParams,'Filename',filenameFig);
end
if pdfs
    for iParams = 1:length(listOfParams)
        makeDelayPlots('pdfs','SelectionParameters',listOfParams(iParams));
    end
end

%% Clean stuff
if cleanDelayPlots
    for iParams = 1:length(listOfParams)
        makeDelayPlots('cleanplots', ...
            'SelectionParameters',listOfParams(iParams));
    end
end
if cleanRadiusDiscAmpl
    listRadiusPlots = dir([figurepath,'*-RadiusDiscAmpl*']);
    cellfun(@(x) delete([figurepath,x]),{listRadiusPlots.name});
end
if cleanPdfs
    for iParams = 1:length(listOfParams)
        makeDelayPlots('cleanpdfs', ...
            'SelectionParameters',listOfParams(iParams));
    end
end

end
