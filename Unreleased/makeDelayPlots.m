function makeDelayPlots(type,varargin)
%MAKEDELAYPLOTS Creates plots from a 3-Channel Series

%% Validate and parse input arguments
p = inputParser;
defaultSelectionParameters = struct('Type','fullcircle', ...
    'Position',[2.5 0.5]);
addParameter(p,'SelectionParameters',defaultSelectionParameters,@isstruct);
parse(p,varargin{:});
c = struct2cell(p.Results);
[selParams] = c{:};

% make all if nothing is specified
if nargin == 0
    type = 'all';
end

%% Find out what needs to be done
[delayMeanVarX,delayDiscAmpl,movieWigner2D,movieWigner3D, ...
    cleanall] = deal(false);
% User request
switch type
    case 'all'
        delayMeanVarX = true;
        delayDiscAmpl = true;
        movieWigner2D = true;
        movieWigner3D = true;
    case 'plots'
        delayMeanVarX = true;
        delayDiscAmpl = true;
    case 'cleanall'
        cleanall = true;
end

% Look what is already there
selStr = selParamsToStr(selParams);
if ~isempty(dir(['*-DelayMeanVarX-',selStr,'*']))
    delayMeanVarX = false;
end
if ~isempty(dir(['*-DelayDiscAmpl-',selStr,'*']))
    delayDiscAmpl = false;
end
if ~isempty(dir(['*-WignerMovie3D-',selStr,'*']))
    movieWigner3D = false;
end
if ~isempty(dir(['*-WignerMovie2D-',selStr,'*']))
    movieWigner2D = false;
end

% Find dependencies that need to be created
[makeTable] = deal(false);
if (delayMeanVarX || delayDiscAmpl) && ...
        isempty(seriesRead3ChTable(selParams))
    makeTable = true;
end

%% Make
datestring = datestr(date,'yyyy-mm-dd');
if makeTable
    T = series3Ch('SelectionParameters',selParams);
else
    T = seriesRead3ChTable(selParams);
end
if delayMeanVarX
    plotSeries3Ch(T,'Type','DelayMeanVarX','Filename', ...
        [datestring,'-DelayMeanVarX-',selStr,'.fig']);
end
if delayMeanVarX
    plotSeries3Ch(T,'Type','DelayDiscAmpl','Filename', ...
        [datestring,'-DelayDiscAmpl-',selStr,'.fig']);
end
if movieWigner2D || movieWigner3D
    series3Ch('SaveWigner',true,'SelectionParameters',selParams);
end
if movieWigner2D
    seriesWignerMovie('Image',true,'Narrow',true);
end
if movieWigner3D
    seriesWignerMovie('Narrow',true);
end


%% Make clean
if cleanall
    listMeanVarX = dir(['*-DelayMeanVarX-',selStr,'*']);
    cellfun(@delete,{listMeanVarX.name});
    listDelayDiscAmpl = dir(['*-DelayDiscAmpl-',selStr,'*']);
    cellfun(@delete,{listDelayDiscAmpl.name});
    listWignerMovie2D = dir(['*-WignerMovie2D-',selStr,'*']);
    cellfun(@delete,{listWignerMovie2D.name});
    listWignerMovie3D = dir(['*-WignerMovie3D-',selStr,'*']);
    cellfun(@delete,{listWignerMovie3D.name});
end

end

