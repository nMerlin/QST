function makeSelectionPlots(type,varargin)
%MAKESELECTIONPLOTS

%% Validate and parse input arguments
p = inputParser;
defaultListOfParams = {};
addParameter(p,'ListOfParams',defaultListOfParams,@isstruct);
defaultPeriod = 2; % number of periods expected in one piezo segment. Important for phase computation. 
addParameter(p,'Period',defaultPeriod,@isnumeric);
defaultRecomputeTheta = false;
addParameter(p,'RecomputeTheta',defaultRecomputeTheta,@islogical);
defaultSavePostselection = false;
addParameter(p,'SavePostselection',defaultSavePostselection,@islogical);
defaultSaveTheta = false;
addParameter(p,'SaveTheta',defaultSaveTheta,@islogical);
defaultRecomputeOrth = false;
addParameter(p,'RecomputeOrth',defaultRecomputeOrth,@islogical);
defaultSaveOrth = false;
addParameter(p,'SaveOrth',defaultSaveOrth,@islogical);
defaultParameter = 'delay';
addParameter(p,'Parameter',defaultParameter,@isstr);
defaultRemoveModulation = false;
addParameter(p,'RemoveModulation',defaultRemoveModulation,@islogical);
defaultRange = 0.3;
addParameter(p,'Range',defaultRange,@isvector);
defaultVaryAPS = false;
addParameter(p,'VaryAPS',defaultVaryAPS,@islogical);
defaultXUnit = 'fs';
addParameter(p,'XUnit',defaultXUnit,@isstr);
defaultFitType = 'gauss';
addParameter(p,'FitType',defaultFitType,@isstr);
defaultChannelAssignment = [1,2,3]; %[target,ps_piezo_fast,ps_piezo_slow]
addParameter(p,'ChannelAssignment',defaultChannelAssignment,@isvector);
defaultCorrRemove = 'yes';
addParameter(p,'CorrRemove',defaultCorrRemove,@isstr);
defaultZeroDelay = 0;
addParameter(p,'ZeroDelay',defaultZeroDelay,@isnumeric);
defaultMeanNs = [13 13 9]; %the mean photon numbers over the total delay series, used for remMod.
% [meanNPsFast,meanNPsSlow,meanNTg]
addParameter(p,'MeanNs',defaultMeanNs,@isvector);
defaultPlotrelative = false; %plots the data relative to the target photon number
addParameter(p,'Plotrelative',defaultPlotrelative,@islogical);
parse(p,varargin{:});
c = struct2cell(p.Results);
[chAssign,corrRemove,fitType,listOfParams,meanNs,parameter,periodsPerSeg,plotrelative,range,recomputeOrth,recomputeTheta,remMod, ...
    saveOrth,saveps,savetheta,varyAPS,xUnit,zeroDelay] = c{:};

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
    cleanRadiusMeanVarSigma,thicknessMeanVarSigma, ...
    cleanThicknessMeanVarSigma,nWithoutPostselection,makeTable] = deal(false);
switch type
    case 'all'
        delayPlots = true;
        radiusDiscAmpl = true;
        radiusMeanVar = true;
        radiusMeanVarSigma = false;
        radiusDiscN = true;
        radiusG2 = true;
        nWithoutPostselection = true;
        pdfs = false;
    case 'table' %only computes the series table without plotting
        makeTable = true;
    case 'G2'
        radiusG2 = true;
    case 'radiusPlots'
        radiusDiscAmpl = true;
        radiusMeanVar = true;
        %radiusMeanVarSigma = true;
        radiusDiscN = true;
        radiusG2 = true;
        nWithoutPostselection = true;
        %delayPlots = true;
    case 'cleanRadiusPlots'
        cleanRadiusDiscAmpl = true;
        cleanRadiusMeanVar = true;
        cleanRadiusMeanVarSigma = true;
        cleanRadiusDiscN = true;
        cleanRadiusG2 = true;
        delayPlots = true;
    case 'thicknessPlots'
        thicknessMeanVar = true;
        thicknessMeanVarSigma = true;
        thicknessMeanVarMin = true;
        thicknessG2 = true;
        delayPlots = true;
    case 'cleanThicknessPlots'
        cleanThicknessMeanVar = true;
        cleanThicknessMeanVarSigma = true;
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
if makeTable
    for iParams = 1:length(listOfParams)
        makeDelayPlots('table','SelectionParameters',listOfParams(iParams),...
        'Period',periodsPerSeg,'RecomputeTheta',recomputeTheta,'RecomputeOrth',recomputeOrth,...
        'SaveOrth',saveOrth,'SavePostselection',saveps,'SaveTheta',savetheta,'Parameter',parameter,...
        'RemoveModulation',remMod,'Range',range,'XUnit',xUnit,'VaryAPS',varyAPS,'ChannelAssignment',...
        chAssign,'CorrRemove',corrRemove,'ZeroDelay',zeroDelay,'MeanNs',meanNs);  
    end
end
if delayPlots
    for iParams = 1:length(listOfParams)
        makeDelayPlots('plots','SelectionParameters',listOfParams(iParams),...
        'Period',periodsPerSeg,'RecomputeTheta',recomputeTheta,'RecomputeOrth',recomputeOrth,...
        'SaveOrth',saveOrth,'SavePostselection',saveps,'SaveTheta',savetheta,'Parameter',parameter,...
        'RemoveModulation',remMod,'Range',range,'XUnit',xUnit,'VaryAPS',varyAPS,'fitType',fitType,...
        'ChannelAssignment',chAssign,'CorrRemove',corrRemove,'ZeroDelay',zeroDelay,'MeanNs',meanNs);  
    end
end
if radiusDiscAmpl
    filenameFig = [figurepath,datestring,'-RadiusDiscAmpl.fig'];
    plotSeriesPostselections(listOfParams,'Filename',filenameFig, ...
        'Type','Amplitude','XUnit',xUnit,'VaryAPS',varyAPS,'fitType',fitType,'RemoveModulation',remMod,'Range',range,'ZeroDelay',zeroDelay,'Plotrelative',plotrelative);
    plotSeriesPostselections(listOfParams,'Filename',[figurepath,datestring,'-RadiusDiscAmplLog.fig'], ...
        'Type','AmplitudeLog','XUnit',xUnit,'VaryAPS',varyAPS,'fitType',fitType,'RemoveModulation',remMod,'Range',range,'ZeroDelay',zeroDelay,'Plotrelative',plotrelative);
end
if radiusMeanVar
    filenameFig = [figurepath,datestring,'-RadiusMeanVar.fig'];
    plotSeriesPostselections(listOfParams,'Filename',filenameFig, ...
        'Type','MeanVar','XUnit',xUnit,'VaryAPS',varyAPS,'fitType',fitType,'RemoveModulation',remMod,'Range',range,'ZeroDelay',zeroDelay,'Plotrelative',plotrelative);
     plotSeriesPostselections(listOfParams,'Filename',[figurepath,datestring,'-VarQ.fig'], ...
        'Type','VarQ','XUnit',xUnit,'VaryAPS',varyAPS,'fitType',fitType,'RemoveModulation',remMod,'Range',range,'ZeroDelay',zeroDelay,'Plotrelative',plotrelative);
     plotSeriesPostselections(listOfParams,'Filename',[figurepath,datestring,'-VarP.fig'], ...
        'Type','VarP','XUnit',xUnit,'VaryAPS',varyAPS,'fitType',fitType,'RemoveModulation',remMod,'Range',range,'ZeroDelay',zeroDelay,'Plotrelative',plotrelative);
end
if radiusMeanVarSigma
    filenameFig = [figurepath,datestring,'-RadiusMeanVarSigma.fig'];
    plotSeriesPostselections(listOfParams,'Filename',filenameFig, ...
        'Type','MeanVarSigma','XUnit',xUnit,'VaryAPS',varyAPS,'fitType',fitType,'RemoveModulation',remMod,'Range',range,'ZeroDelay',zeroDelay,'Plotrelative',plotrelative);   
end
if radiusDiscN
    filenameFig = [figurepath,datestring,'-RadiusDiscN.fig'];
    plotSeriesPostselections(listOfParams,'Filename',filenameFig, ...
        'Type','DiscN','XUnit',xUnit,'VaryAPS',varyAPS,'fitType',fitType,'RemoveModulation',remMod,'Range',range,'ZeroDelay',zeroDelay,'Plotrelative',plotrelative);
end
if radiusG2
    filenameFig = [figurepath,datestring,'-RadiusG2.fig'];
    plotSeriesPostselections(listOfParams,'Filename',filenameFig, ...
        'Type','G2','XUnit',xUnit,'VaryAPS',varyAPS,'fitType',fitType,'RemoveModulation',remMod,'Range',range,'ZeroDelay',zeroDelay,'Plotrelative',plotrelative);
end    
if nWithoutPostselection
    filenameFig = [figurepath,datestring,'-nWithoutPostselection.fig'];
    plotSeriesPostselections(listOfParams,'Filename',filenameFig, ...
        'Type','nWithoutPostselection','XUnit',xUnit,'VaryAPS',varyAPS,'fitType',fitType,'RemoveModulation',remMod,'Range',range,'ZeroDelay',zeroDelay,'Plotrelative',plotrelative);
end   
if thicknessMeanVar
    filenameFig = [figurepath,datestring,'-ThicknessMeanVar.fig'];
    plotSeriesPostselections(listOfParams,'Filename',filenameFig, ...
        'Type','ThicknessMeanVar','XUnit',xUnit,'VaryAPS',varyAPS,'fitType',fitType,'RemoveModulation',remMod,'Range',range,'ZeroDelay',zeroDelay,'Plotrelative',plotrelative);
end
if thicknessMeanVarSigma
    filenameFig = [figurepath,datestring,'-ThicknessMeanVarSigma.fig'];
    plotSeriesPostselections(listOfParams,'Filename',filenameFig, ...
        'Type','ThicknessMeanVarSigma','XUnit',xUnit,'VaryAPS',varyAPS,'fitType',fitType,'RemoveModulation',remMod,'Range',range,'ZeroDelay',zeroDelay,'Plotrelative',plotrelative);
end
if thicknessMeanVarMin
    filenameFig = [figurepath,datestring,'-ThicknessMeanVarMin.fig'];
    plotSeriesPostselections(listOfParams,'Filename',filenameFig, ...
        'Type','ThicknessMeanVarMin','XUnit',xUnit,'VaryAPS',varyAPS,'fitType',fitType,'RemoveModulation',remMod,'Range',range,'ZeroDelay',zeroDelay,'Plotrelative',plotrelative);
end
if thicknessG2
    filenameFig = [figurepath,datestring,'-ThicknessG2.fig'];
    plotSeriesPostselections(listOfParams,'Filename',filenameFig, ...
        'Type','ThicknessG2','XUnit',xUnit,'VaryAPS',varyAPS,'fitType',fitType,'RemoveModulation',remMod,'Range',range,'ZeroDelay',zeroDelay,'Plotrelative',plotrelative);
end
if pdfs
    for iParams = 1:length(listOfParams)
        makeDelayPlots('pdfs','SelectionParameters',listOfParams(iParams),'XUnit',xUnit,'VaryAPS',varyAPS,'RemoveModulation',remMod,'Range',range,'ZeroDelay',zeroDelay,'Plotrelative',plotrelative);
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
if cleanThicknessMeanVarSigma
    list = dir([figurepath,'*-ThicknessMeanVarSigma*']);
    cellfun(@(x) delete([figurepath,x]),{list.name});
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
