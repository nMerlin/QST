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
parse(p,varargin{:});
c = struct2cell(p.Results);
[image,narrow] = c{:};

%% Check data availability and capture data table
dispstat('','init')
dispstat('seriesWignerMovie: Checking folder ''post-data''', ...
    'timestamp','keepthis');
T = seriesRead3ChTable();
if ~isempty(T) && ismember('Delay',T.Properties.VariableNames)
    files = T.Filename;
    delays = T.Delay;
    titles = cell(length(delays),1);
    for iDelay = 1:length(delays)
        titles{iDelay} = ['Delay: ',num2str(delays(iDelay)),' fs'];
    end
else
    titles = {};
    %% Discover *.mat files
    if exist('post-data','dir')
        filestruct = dir('post-data/*.mat');
        files = {filestruct.name};
    else
        warning('Could not find the ''post-data'' directory!');
    end
end

%% Capture the available WF variables
allWF = {};
for iFile = 1:length(files)
    ss = strsplit(files{iFile},'.');
    dispstat(['seriesWignerMovie: Loading ','post-data/',ss{1}, ...
        '-postselection.mat'],'timestamp');
    load(['post-data/',ss{1},'-postselection.mat']);
    if exist('WF','var')
        allWF(end+1) = {WF};
        clear WF;
    else
        titles{iFile} = {};
    end
end
titles = titles(~cellfun(@isempty,titles));

%% Generate Movie
plotfun = @(wigfun,params) plotWigner(wigfun,params,varargin{:});
plotMovie(plotfun,allWF,'Titles',titles);

end
