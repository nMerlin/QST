function plotMovie(plotFun,plotInputs,varargin)
%PLOTMOVIE Uses PLOTFUN to make a movie of PLOTINPUTS
%   Detailed explanation goes here

%% Validate and parse input arguments
p = inputParser;
defaultDelays = 24;
addParameter(p,'Delays',defaultDelays,@isnumeric);
defaultFigureVisible = 'off';
addParameter(p,'FigureVisible',defaultFigureVisible,@isstr);
defaultFilename = 'Movie.mp4';
addParameter(p,'Filename',defaultFilename,@isstr);
defaultTitles = {};
addParameter(p,'Titles',defaultTitles,@iscell);
defaultZLim = [NaN NaN];
addParameter(p,'ZLim',defaultZLim,@isvector);
parse(p,varargin{:});
c = struct2cell(p.Results);
[delays,figurevisible,filename,titles,zlim] = c{:};

%% Preparing figure
dispstat('','init')
dispstat(['Working on ',filename,' ...'],'timestamp','keepthis');
if ~isempty(findobj('type','figure'))
    delete(findall(0)); % close all figures
end
%fig = figure('Visible',figurevisible);
fig = figure;
axis tight;  % set axis limit to the range of the data
axis manual;
movie = VideoWriter(filename, 'MPEG-4');
movie.Quality = 100;
open(movie);

plotFunParams.Handle = get(fig,'CurrentAxes');
for iInput = 1:length(plotInputs)
    dispstat(['Generating frame ',num2str(iInput),' of ', ...
        num2str(length(plotInputs)),'.'],'timestamp');
    delete(findall(gcf,'Type','light'));
    plotFun(plotInputs{iInput},plotFunParams);
    if ~isnan(zlim(1))
        set(gca,'ZLim',zlim);
    end
    if isempty(titles)
        title(['Frame ',num2str(iInput),' of ', ...
            num2str(length(plotInputs))]);
    else
        title(titles{iInput});
    end
    frame = print(fig,'-r150','-RGBImage');
    for iDelay = 1:delays
        writeVideo(movie, frame);
    end % iDelay
end % iInput

end

