function seriesWignerMovie()
%SERIESWIGNERMOVIE Creates a movie from series of Wigner functions. Is a
%wrapper function for 'plotMovie'.
%
% Prerequisites:
%   - expects 'WF' variable in the *.mat files in a 'post-data' folder

%% Discover *.mat files
dispstat('','init')
dispstat('seriesWignerMovie: Checking folder ''post-data''', ...
    'timestamp','keepthis');
if exist('post-data','dir')
    filestruct = dir('post-data/*.mat');
    files = {filestruct.name};
else
    warning('Could not find the ''post-data'' directory!');
end

%% Capture the available WF variables
allWF = {};
for iFile = 1:length(files)
    dispstat(['seriesWignerMovie: Loading ','post-data/',files{iFile}], ...
        'timestamp');
    load(['post-data/',files{iFile}]);
    if exist('WF','var')
        allWF(end+1) = {WF};
        clear WF;
    end
end

%% Generate Movie
plotMovie(@(wigfun,params) plotWigner(wigfun,params),allWF);

end
