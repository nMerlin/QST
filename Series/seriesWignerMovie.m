function seriesWignerMovie(varargin)
%SERIESWIGNERMOVIE Creates a movie from series of Wigner functions. Is a
%wrapper function for 'plotMovie'.
%
% Optional Input Arguments:
%   'Image': Default is 'false' (3D plot). If 'true',
%       create a 2D plot of WF.
%   'Narrow': Default is 'false'. Narrower axis limits.
%
% Prerequisites:
%   - expects 'WF' variable in the *.mat files in a 'post-data' folder

%% Validate and parse input arguments
p = inputParser;
defaultImage = false;
addParameter(p,'Image',defaultImage,@islogical);
defaultNarrow = false;
addParameter(p,'Narrow',defaultNarrow,@islogical);
defaultSelParams = struct('Type','fullcircle','Position',[2.5,0.5]);
addParameter(p,'SelectionParameters',defaultSelParams,@isstruct);
parse(p,varargin{:});
c = struct2cell(p.Results);
[image,narrow,selParams] = c{:};
wigParams.Image = image;
wigParams.Narrow = narrow;

%% Check data availability and capture data table
dispstat('','init')
dispstat('seriesWignerMovie: Checking folder ''post-data''', ...
    'timestamp','keepthis');
T = seriesRead3ChTable();
if ~isempty(T) && ismember('Delay',T.Properties.VariableNames)
    files = T.Filename;
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
    ss = strsplit(files{iFile},'.');
    dispstat(['seriesWignerMovie: Loading ','post-data/',ss{1}, ...
        '-',selStr,'.mat'],'timestamp');
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

%% Generate Movie
plotfun = @(wigfun,params) plotWigner(wigfun,params,wigParams);

filename = [datestr(date,'yyyy-mm-dd-'),'WignerMovie', ...
    num2str(2+not(image)),'D-',selStr,'.mp4'];
plotMovie(plotfun,allWF,'Titles',titles,'Filename',filename);

end
