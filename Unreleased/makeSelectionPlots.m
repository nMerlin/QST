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
    cleanRadiusDiscN,thicknessMeanVar,cleanThicknessMeanVar, ...
    thicknessMeanVarMin,cleanThicknessMeanVarMin,radiusG2,cleanRadiusG2,...
    thicknessG2,cleanThicknessG2,radiusMeanVarSigma, ...
    cleanRadiusMeanVarSigma] = deal(false);
switch type
    case 'all'
        delayPlots = true;
        radiusDiscAmpl = true;
        radiusMeanVar = true;
        radiusMeanVarSigma = true;
        radiusDiscN = true;
        radiusG2 = true;
        pdfs = true;
    case 'G2'
        radiusG2 = true;
    case 'radiusPlots'
        radiusDiscAmpl = true;
        radiusMeanVar = true;
        radiusMeanVarSigma = true;
        radiusDiscN = true;
        radiusG2 = true;
        delayPlots = true;
    case 'cleanRadiusPlots'
        cleanRadiusDiscAmpl = true;
        cleanRadiusMeanVar = true;
        cleanRadiusMeanVarSigma = true;
        cleanRadiusDiscN = true;
        cleanRadiusG2 = true;
        delayPlots = true;
    case 'thicknessPlots'
        thicknessMeanVar = true;
        thicknessMeanVarMin = true;
        thicknessG2 = true;
        delayPlots = true;
    case 'cleanThicknessPlots'
        cleanThicknessMeanVar = true;
        cleanThicknessMeanVarMin = true;
        cleanThicknessG2 = true;
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
if radiusMeanVarSigma
    filenameFig = [figurepath,datestring,'-RadiusMeanVarSigma.fig'];
    plotSeriesPostselections(listOfParams,'Filename',filenameFig, ...
        'Type','MeanVarSigma');
end
if radiusDiscN
    filenameFig = [figurepath,datestring,'-RadiusDiscN.fig'];
    plotSeriesPostselections(listOfParams,'Filename',filenameFig, ...
        'Type','DiscN');
end
if radiusG2
    filenameFig = [figurepath,datestring,'-RadiusG2.fig'];
    plotSeriesPostselections(listOfParams,'Filename',filenameFig, ...
        'Type','G2');
end
if thicknessMeanVar
    filenameFig = [figurepath,datestring,'-ThicknessMeanVar.fig'];
    plotSeriesPostselections(listOfParams,'Filename',filenameFig, ...
        'Type','ThicknessMeanVar');
end
if thicknessMeanVarMin
    filenameFig = [figurepath,datestring,'-ThicknessMeanVarMin.fig'];
    plotSeriesPostselections(listOfParams,'Filename',filenameFig, ...
        'Type','ThicknessMeanVarMin');
end
if thicknessG2
    filenameFig = [figurepath,datestring,'-ThicknessG2.fig'];
    plotSeriesPostselections(listOfParams,'Filename',filenameFig, ...
        'Type','ThicknessG2');
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
if cleanRadiusMeanVarSigma
    listRadiusPlots = dir([figurepath,'*-RadiusMeanVarSigma*']);
    cellfun(@(x) delete([figurepath,x]),{listRadiusPlots.name});
end
if cleanRadiusDiscN
    listRadiusPlots = dir([figurepath,'*-RadiusDiscN*']);
    cellfun(@(x) delete([figurepath,x]),{listRadiusPlots.name});
end
if cleanRadiusG2
    listRadiusG2 = dir([figurepath,'*-RadiusG2*']);
    cellfun(@(x) delete([figurepath,x]),{listRadiusG2.name});
end
if cleanThicknessMeanVar
    listThicknessMeanVar = dir([figurepath,'*-ThicknessMeanVar*']);
    cellfun(@(x) delete([figurepath,x]),{listThicknessMeanVar.name});
end
if cleanThicknessMeanVarMin
    listThicknessMeanVarMin = dir([figurepath,'*-ThicknessMeanVarMin*']);
    cellfun(@(x) delete([figurepath,x]),{listThicknessMeanVarMin.name});
end
if cleanThicknessG2
    listThicknessG2 = dir([figurepath,'*-ThicknessG2*']);
    cellfun(@(x) delete([figurepath,x]),{listThicknessG2.name});
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
