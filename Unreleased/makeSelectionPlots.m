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

[radiusDiscAmpl] = deal(false);
switch type
    case 'all'
        radiusDiscAmpl = true;
end

%% Resolve dependencies

%% Make stuff
if radiusDiscAmpl
    plotSeriesPostselections(listOfParams);
end

end
