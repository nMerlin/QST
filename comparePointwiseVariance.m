function [] = comparePointwiseVariance( varargin )
%This function plots the pointwise variance and locs of a detector standard test.
%   PLOTSHOTNOISE('verbose'): Shows log output.
%   The script assumes that the datafiles are in the 'raw-data' directory.
%   The filename-convention is '03-0.1mW-*.raw'.
%
%   See also LOAD8BITBINARY.

% Optional input arguments
verbose = 0;
quiet = 'notquiet';
if nargin > 0
    for i = 1:nargin
        eval([varargin{i} '=1;']);
    end
end
if verbose == 0
    quiet = 'quiet';
end

% Parameters & Variables
outputFilename1 = 'pointwise-variance-plot-all.jpg';
%outputFilename1 = 'pointwise-variance-plot-first3.jpg';
outputFilename2 = 'Nlocs.jpg';
outputFiletype = '-djpeg';
range = 275; %index range for pointwise variance plot
dataStruct = struct('filename',{},'powerLO',{}, 'Nlocs', {});

%%% Create data overview
dispstat('','init',quiet);
dispstat('Checking filenames ...','timestamp','keepthis',quiet);
rawDataContents = dir('raw-data');
for name = {rawDataContents.name}
    % Loop only over *.raw files
    filename = cell2mat(name);
    if not(isempty(regexpi(filename,'.raw.','match'))) || isempty(regexpi(filename,'.raw','match'))
        continue
    end
    
    % Fetch number of measurement
    numberToken = regexpi(filename,'^([0123456789]*)-','tokens');
    number = str2double(cell2mat(numberToken{1}));
    dataStruct(number).filename = filename;
    
    % Get LO power
    powerToken = regexpi(filename,'-([0123456789.]*)mW','tokens');
    dataStruct(number).powerLO = str2double(cell2mat(powerToken{1}));
end


% Calculate and plot pointwise variance for each LO-power
dispstat('Calculating and plotting ...','timestamp','keepthis',quiet);
close all;
figure;

maxvar0 = 0;

for number=1:size(dataStruct,2)
%for number = 1:3
    [data8bit,config,~] = load8BitBinary(dataStruct(number).filename,'dontsave');
    [locs,pvar] = pointwiseVariance(data8bit);
    maxvar = max(pvar);
    if maxvar > maxvar0
        maxvar0 = maxvar;
    end
    dataStruct(number).Nlocs = size(locs,1);
    time = 1/config.SpectrumCard.Clock.SamplingRate0x28MHz0x29_DBL;
    t = (1:range) *time;
    fullvar = var(double(data8bit(:)));
    plot(t, pvar(1:range),'DisplayName', strcat(num2str(dataStruct(number).powerLO), ' mW '));
    %plot(t, pvar(1:range),'Color',[1-1/number 1/number 1/number],'DisplayName',...
     %   strcat(num2str(dataStruct(number).powerLO),' mW ; Var = ', num2str(fullvar)));  %writes the variance of the full data in the legend
    hold on;
    for i = 1:find(locs>=range,1)
        p = plot(locs(i)*time,pvar(locs(i)),'o','MarkerFaceColor', 'k', 'MarkerEdgeColor','k', ...
    'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'locs');
        set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  %excludes this from legend
    end
    hold on;
end
%hold off;
legend( 'Location','northeast');
axis ([0 range*time 0 maxvar0 + 1]);
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultTextInterpreter','latex');
xlabel('Time [s]');
ylabel('Pointwise Variance [a. u.]')

% Saving figure
print(outputFilename1,outputFiletype);

%plot Number of locs over Power

hold off;
close all;
figure;
powerLO = cell2mat({dataStruct.powerLO});
Nlocs = cell2mat({dataStruct.Nlocs});
disp(powerLO);
disp(Nlocs);
plot(powerLO, Nlocs,'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor','k', ...
    'LineWidth', 2, 'MarkerSize', 4);
xlabel('P_{LO} [mW]');
ylabel('Number of Locs')
axis ([0 max(powerLO)+1 min(Nlocs)-2 max(Nlocs)+2]);
print(outputFilename2,outputFiletype);


dispstat('Finished!','timestamp',quiet);

end