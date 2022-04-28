function [] = SwitchingPhotonnumberSeries( level_on, level_off, varargin )
%% Validate and parse input arguments
p = inputParser;
defaultSmooth = false; % whether data is smoothed
addParameter(p,'Smooth',defaultSmooth,@islogical);
defaultFolder = 'g2data-X1-Concatenated'; % which folder contain the data
addParameter(p,'Folder',defaultFolder);
parse(p,varargin{:});
c = struct2cell(p.Results);
[folder, smoothing] = c{:};
%%
targetfolder = [folder '-SwitchingPlots'];
if ~exist(targetfolder,'dir')
    mkdir(targetfolder)
end
filenameOptions = ['-smooth-' num2str(smoothing) '-level_on-' num2str(level_on) '-level_off-' num2str(level_off)];

[filenames,numbers,Is]= getParametersFromFilenames('Folder',folder);
[period_mean_ons,period_mean_offs,fliprates] = deal(zeros(length(filenames),1));

for fileI = 1:length(filenames)
    load([folder '\' cell2mat(filenames(fileI))], 'times', 'ada', 'g2vec');
    plotFilename = [targetfolder '\' cell2mat(filenames(fileI)) ];
    [ ~,~, period_mean_on, period_mean_off, fliprate ] = ...
        SwitchingPhotonnumber( times, ada, g2vec, level_on, level_off, plotFilename,'Smooth', smoothing );
    period_mean_ons(fileI) = period_mean_on;
    period_mean_offs(fileI) = period_mean_off;
    fliprates(fileI) = fliprate ;
end

[~,units] = convenientUnits(times,'s');
save([targetfolder '\SwitchingResults' filenameOptions '.mat'],...
    'filenames','numbers','Is','period_mean_ons','period_mean_offs','fliprates');

%% plot fliprates 
plot(Is,fliprates,'o-');
graphicsSettings;
xlabel('Power (mW)');
ylabel('Flip rate (Hz)');
savefig([targetfolder '\Fliprates' filenameOptions '.fig']);
print([targetfolder '\Fliprates' filenameOptions '.png'],'-dpng');
close;

%% plot mean dwell times
[AX,H1,H2] = plotyy(Is,period_mean_ons,Is,period_mean_offs);
H1.Marker = 'o';
H2.Marker = '*';
xlabel('Power (mW)');
ylabel(['Time (' units ')']);
ylabel(AX(2),['Time (' units ')']);
legend('"On" times','"Off" times','Location','best');
set(H1,'linewidth',2);
    set(AX(1),'linewidth',2,...
        'FontSize',20,'FontName','Arial',...
        'TickDir','in');
set(H2,'linewidth',2);
    set(AX(2),'linewidth',2,...
        'FontSize',20,'FontName','Arial',...
        'TickDir','in');
    fig = figure(1);              
set(fig,'Color','w','Units','centimeters','Position',[1,1,45,30],'PaperPositionMode','auto');
savefig([targetfolder '\MeanTimesVsPower' filenameOptions '.fig']);
print([targetfolder '\MeanTimesVsPower' filenameOptions '.png'],'-dpng');
close;

% plot ratio
plot(Is,period_mean_ons./period_mean_offs,'-o');
xlabel('Power (mW)');
ylabel('t_{on} / t_{off}');
graphicsSettings;
savefig([targetfolder '\RatioTimesVsPower' filenameOptions '.fig']);
print([targetfolder '\RatioTimesVsPower' filenameOptions '.png'],'-dpng');
close;

end