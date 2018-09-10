function seriesWignerMovie(varargin)
%SERIESWIGNERMOVIE Creates a movie from series of Wigner functions. Is a
%wrapper function for 'plotMovie'.
%
% Optional Input Arguments:
%   'Dimensions': Specify '2D' or '3D' movie. Default is 'both'.
%   'Image': Default is 'false' (3D plot). If 'true',
%       create a 2D plot of WF.
%   'Narrow': Default is 'false'. Narrower axis limits.
%
% Prerequisites:
%   - expects 'WF' variable in the *.mat files in a 'post-data' folder

%% Validate and parse input arguments
p = inputParser;
defaultDimensions = 'both';
addParameter(p,'Dimensions',defaultDimensions,@isstr);
defaultPQ = -20:0.125:20;
addParameter(p,'PQ',defaultPQ,@isvector);
defaultSelParams = struct('Type','fullcircle','Position',[2.5,0.5]);
addParameter(p,'SelectionParameters',defaultSelParams,@isstruct);
defaultZLim = [NaN NaN];
addParameter(p,'ZLim',defaultZLim,@isvector);
parse(p,varargin{:});
c = struct2cell(p.Results);
[dimensions,pq,selParams,zlimits] = c{:};
wigParams.PQ = pq;
wigParams.ZLim = zlimits;

%% Check data availability and capture data table
dispstat('','init')
dispstat('Creating movies ...','timestamp','keepthis');
T = seriesRead3ChTable();
if ~isempty(T) && ismember('Delay',T.Properties.VariableNames)
    files = T.Filename;
    % Remove filetype
    [~,files] = cellfun(@fileparts,files,'UniformOutput',false);
    delays = T.Delay;
else
    %% Discover *.mat files
    if exist('post-data','dir')
        filestruct = dir('post-data/*.mat');
        files = {filestruct.name};
    else
        warning('Could not find the ''post-data'' directory!');
    end
end

%% Capture the available WF variables
selStr = selParamsToStr(selParams);
allWF = {};
for iFile = 1:length(files)
    ss = strsplit(files{iFile},'.mat');
    dispstat(['Loading ','post-data/', ...
        ss{1},'-',selStr,'.mat'],'timestamp');
    load(['post-data/',ss{1},'-',selStr,'.mat']);
    if exist('WF','var')
        allWF(end+1) = {WF};
        clear WF;
    else
        delays(iFile) = NaN;
    end
end
delays = delays(~isnan(delays));

%% Sort for delays
[delays,iDelays] = sort(delays);
allWF = allWF(iDelays);

%% Create meaningful plot titles
if ~isempty(delays)
    titles = cell(length(delays),1);
    for iDelay = 1:length(delays)
        titles{iDelay} = ['Delay: ',num2str(round(delays(iDelay))),' fs'];
    end
else
    titles = {};
end

%% Generate Movie(s)
[twoD,threeD] = deal(false);
if strcmp(dimensions,'both')
    twoD = true;
    threeD = true;
elseif strcmp(dimensions,'2D')
    twoD = true;
elseif strcmp(dimensions,'3D')
    threeD = true;
end

if twoD
    wigParams.Style = '2D';
    plotfun = @(wigfun,params) plotWigner(wigfun,params,wigParams);
    filename = [datestr(date,'yyyy-mm-dd-'),'WignerMovie-', ...
        wigParams.Style,'-',selStr,'.mp4'];
    plotMovie(plotfun,allWF,'Titles',titles,'Filename',filename);
end
if threeD
    % compute colormap
    minVal = min(cell2mat(cellfun(@(A) min(min(A)),allWF, ...
        'UniformOutput',false)));
    maxVal = max(cell2mat(cellfun(@(A) max(max(A)),allWF, ...
        'UniformOutput',false)));
    cmap = hsv(128);
    shift = ceil(abs(minVal)/(maxVal-minVal)*length(cmap));
    cmap = circshift(cmap,-shift,1);
    wigParams.Style = 'advanced';
    wigParams.EdgeColor = 'None';
    wigParams.ColorMap = cmap;
    wigParams.ColorLimits = [minVal,maxVal];
    plotfun = @(wigfun,params) plotWigner(wigfun,params,wigParams);
    filename = [datestr(date,'yyyy-mm-dd-'),'WignerMovie-', ...
        wigParams.Style,'-',selStr,'.mp4'];
    plotMovie(plotfun,allWF,'Titles',titles,'Filename',filename);
end


end
