function makeSelectionPlots(type,varargin)
%MAKESELECTIONPLOTS

%% Validate and parse input arguments
p = inputParser;
defaultListOfParams = {};
addParameter(p,'ListOfParams',defaultListOfParams,@isstruct);
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

[radiusDiscAmpl,cleanRadiusDiscAmpl,cleanDelayPlots,delayPlots,pdfs, ...
    cleanpdfs,radiusMeanVar,cleanRadiusMeanVar,radiusDiscN, ...
    cleanRadiusDiscN,thicknessMeanVar,cleanThicknessMeanVar] = deal(false);
switch type
    case 'all'
        delayPlots = true;
        radiusDiscAmpl = true;
        radiusMeanVar = true;
        radiusDiscN = true;
        pdfs = true;
    case 'radiusPlots'
        radiusDiscAmpl = true;
        radiusMeanVar = true;
        radiusDiscN = true;
        delayPlots = true;
    case 'cleanRadiusPlots'
        cleanRadiusDiscAmpl = true;
        cleanRadiusMeanVar = true;
        cleanRadiusDiscN = true;
        delayPlots = true;
    case 'thicknessPlots'
        thicknessMeanVar = true;
        delayPlots = true;
    case 'cleanThicknessPlots'
        cleanThicknessMeanVar = true;
        delayPlots = true;
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
    plotSeriesPostselections(listOfParams,'Filename',filenameFig, ...
        'Type','Amplitude');
end
if radiusMeanVar
    filenameFig = [figurepath,datestring,'-RadiusMeanVar.fig'];
    plotSeriesPostselections(listOfParams,'Filename',filenameFig, ...
        'Type','MeanVar');
end
if radiusDiscN
    filenameFig = [figurepath,datestring,'-RadiusDiscN.fig'];
    plotSeriesPostselections(listOfParams,'Filename',filenameFig, ...
        'Type','DiscN');
end
if thicknessMeanVar
    filenameFig = [figurepath,datestring,'-ThicknessMeanVar.fig'];
    plotSeriesPostselections(listOfParams,'Filename',filenameFig, ...
        'Type','ThicknessMeanVar');
end
if pdfs
    for iParams = 1:length(listOfParams)
        makeDelayPlots('pdfs','SelectionParameters',listOfParams(iParams));
    end
    % Radius pdfs
    listOfFigures = dir([figurepath,'*-Radius*.fig']);
    [~,figNames] = cellfun(@fileparts,{listOfFigures.name}, ...
        'UniformOutput',false);
    listOfPdfs = dir([pdfpath,'*-Radius*.pdf']);
    [~,pdfNames] = cellfun(@fileparts,{listOfPdfs.name}, ...
        'UniformOutput',false);
    cellfun(@(x) makePdf([figurepath,x,'.fig'],pdfpath), ...
        setdiff(figNames,pdfNames));
    % Thickness pdfs
    listOfFigures = dir([figurepath,'*-Thickness*.fig']);
    [~,figNames] = cellfun(@fileparts,{listOfFigures.name}, ...
        'UniformOutput',false);
    listOfPdfs = dir([pdfpath,'*-Thickness*.pdf']);
    [~,pdfNames] = cellfun(@fileparts,{listOfPdfs.name}, ...
        'UniformOutput',false);
    cellfun(@(x) makePdf([figurepath,x,'.fig'],pdfpath), ...
        setdiff(figNames,pdfNames));
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
if cleanRadiusMeanVar
    listRadiusPlots = dir([figurepath,'*-RadiusMeanVar*']);
    cellfun(@(x) delete([figurepath,x]),{listRadiusPlots.name});
end
if cleanRadiusDiscN
    listRadiusPlots = dir([figurepath,'*-RadiusDiscN*']);
    cellfun(@(x) delete([figurepath,x]),{listRadiusPlots.name});
end
if cleanThicknessMeanVar
    listThicknessMeanVar = dir([figurepath,'*-ThicknessMeanVar*']);
    cellfun(@(x) delete([figurepath,x]),{listThicknessMeanVar.name});
end
if cleanpdfs
    for iParams = 1:length(listOfParams)
        makeDelayPlots('cleanpdfs', ...
            'SelectionParameters',listOfParams(iParams));
    end
    listRadiusPdfs = dir([pdfpath,'*-Radius*']);
    cellfun(@(x) delete([pdfpath,x]),{listRadiusPdfs.name});
    listRadiusPdfs = dir([pdfpath,'*-Thickness*']);
    cellfun(@(x) delete([pdfpath,x]),{listRadiusPdfs.name});
end

end
