function makeSelectionPlots(type,varargin)
%MAKESELECTIONPLOTS

%% Validate and parse input arguments
p = inputParser;
defaultListOfParams = '';
addParameter(p,'ListOfParams',defaultListOfParams,@isstr);
parse(p,varargin{:});
c = struct2cell(p.Results);
[listOfParams] = c{:};

if nargin == 0
    type = 'all';
end
if isempty(listOfParams)
    listMeanVarX = dir('*-series3Ch-*');
    listOfParams = cellfun(@selStrToParams,{listMeanVarX.name});
end

[radiusDiscAmpl,cleanDelayPlots,delayPlots] = deal(false);
switch type
    case 'all'
        radiusDiscAmpl = true;
    case 'cleanDelayPlots'
        cleanDelayPlots = true;
    case 'delayPlots'
        delayPlots = true;
end

%% Resolve dependencies

%% Make stuff
if radiusDiscAmpl
    plotSeriesPostselections(listOfParams);
end
if delayPlots
    for iParams = 1:length(listOfParams)
        makeDelayPlots('plots','SelectionParameters',listOfParams(iParams));
    end
end

%% Clean stuff
if cleanDelayPlots
    for iParams = 1:length(listOfParams)
        makeDelayPlots('cleanplots', ...
            'SelectionParameters',listOfParams(iParams));
    end
end

end
