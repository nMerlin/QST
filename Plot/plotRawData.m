function [] = plotRawData( varargin )
%This function plots the raw data of a detector standard test.
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
outputFilename1 = 'Raw-data.jpg';
outputFiletype = '-djpeg';
range = 2000; %index range for plot
ymax = 10;
ymin = -40;
dataStruct = struct('filename',{},'powerLO',{});

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


% plot 
dispstat('plotting ...','timestamp','keepthis',quiet);
close all;
figure;

for number=1:size(dataStruct,2)
    [data8bit,config,~] = load8BitBinary(dataStruct(number).filename,'dontsave');
    time = 1/config.SpectrumCard.Clock.SamplingRate0x28MHz0x29_DBL * 10^-6 *10^9;  %time in nanoseconds
    t = (1:range) *time;
    plot(t, data8bit(1:range),'DisplayName', strcat(num2str(dataStruct(number).powerLO), ' mW '));
    hold on;
end
hold off;
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultTextInterpreter','latex');
legend( 'Location','northeast');
axis ([0 range*time ymin ymax]);
xlabel('Time [ns]');
ylabel('Fluctuation [a. u.]')

% Saving figure
print(outputFilename1,outputFiletype);

dispstat('Finished!','timestamp',quiet);

end