function plotMovie(plotFun,plotInputs,varargin)
%PLOTMOVIE Uses PLOTFUN to make a movie of PLOTINPUTS
%   Detailed explanation goes here

%% Validate and parse input arguments
p = inputParser;
defaultFilename = 'Movie.mp4';
defaultDelays = 1;
defaultZLim = [NaN NaN];
addParameter(p,'Filename',defaultFilename,@isstr);
addParameter(p,'Delays',defaultDelays,@isnumeric);
addParameter(p,'ZLim',defaultZLim,@isvector);
parse(p,varargin{:});
c = struct2cell(p.Results);
[delays,filename,zlim] = c{:};

h = figure;
axis tight;  % set axis limit to the range of the data
axis manual;
movie = VideoWriter(filename, 'MPEG-4');
movie.Quality = 100;
open(movie);

for iInput = 1:length(plotInputs)
    plotFun(plotInputs{iInput},'Handle',gca);
    if ~isnan(zlim(1))
        set(gca,'ZLim',zlim);
    end
    title(['Iteration ',num2str(iInput),' of ', ...
        num2str(length(plotInputs))]);
    frame = print(h,'-r150','-RGBImage');
    for iDelay = 1:delays
        writeVideo(movie, frame);
    end % iDelay
end % iInput

end

