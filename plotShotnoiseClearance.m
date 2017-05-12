function [  ] = plotShotnoiseClearance( varargin )
%PLOTSHOTNOISECLEARANCE

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
dataStruct = struct('filename',{},'powerLO',{},'filtered',{},...
    'singlediode',{},'noise',{});
fcount = 0;
outputFilename = 'shotnoise-clearance.jpg';
outputFiletype = '-djpeg';
clearanceAveragingLimits = [20 40]; % in MHz
yMin = -85;
yMax = -65;

%%% Create data overview
dispstat('','init',quiet);
dispstat('Checking filenames ...','timestamp','keepthis',quiet);
rawDataContents = dir('raw-data');
for name = {rawDataContents.name}
    % Loop only over *-FFT.txt files & get filename
    filename = cell2mat(name);
    if isempty(regexpi(filename,'FFT.txt','match'))
        continue
    end
    
    fcount = fcount + 1;
    dataStruct(fcount).filename = filename;
    
    % Get LO power
    powerToken = regexpi(filename,'-([0123456789.]*)mW','tokens');
    dataStruct(fcount).powerLO = str2double(cell2mat(powerToken{1}));
    
    % Pure or filtered?
    dataStruct(fcount).filtered = not(isempty(regexpi(filename,'filtered','tokens')));
    
    % Single diode measurement?
    dataStruct(fcount).singlediode = not(isempty(regexpi(filename,'single','tokens')));
end

dispstat('Loading data for plots ...','timestamp','keepthis',quiet);
%%% Process data for plotting
powerLO = cell2mat({dataStruct.powerLO});
[maxPower,~] = max(powerLO);
[minPower,~] = min(powerLO);
dispstat(['MaxPower: ' num2str(maxPower) '; MinPower: ' num2str(minPower)], ...
    'timestamp','keepthis',quiet);

cd('raw-data');

% Read header information and calculate frequency values
dispstat('Reading header information ...','timestamp','keepthis',quiet);
fileID = fopen(dataStruct(1).filename,'r');
header = textscan(fileID,'%s',11,'delimiter','\n');
fclose(fileID);
header = strjoin(header{1});
headerPost = regexpi(header,'s (\d*) \(Post\)','tokens');
headerPost = str2double(cell2mat(headerPost{1}));
headerSampleRate = regexpi(header,...
    's (\d*) \(SampleRate\)','tokens');
headerSampleRate = str2double(cell2mat(headerSampleRate{1}));
frequencyAxis = 1:headerPost;
frequencyAxis = frequencyAxis.*...
    ((1/(headerSampleRate/10000000.0/headerPost/1000))/2)/headerPost;

dispstat('Calculate shot noise clearance ...','timestamp','keepthis',quiet);
startIndex = find(frequencyAxis > clearanceAveragingLimits(1)*1000000,1);
stopIndex = find(frequencyAxis > clearanceAveragingLimits(2)*1000000,1);
for k=1:size(dataStruct,2)
    data = dlmread(dataStruct(k).filename,' ',12,0);
    data(:,1) = frequencyAxis;
    dataStruct(k).noise = mean(data(startIndex:stopIndex,2));
end
cd('..');

noise = cell2mat({dataStruct.noise});
ft = fittype({'log(x)','1'});
logFit = fit(powerLO(5:end)',noise(5:end)',ft);
fitY = logFit(powerLO);

% Plot noise power for 1mW to 12mW
semilogx(powerLO,noise,'o');
text(powerLO(9),noise(9)-0.5,'5 mW');
axis([minPower maxPower yMin yMax]);
xlabel('LO Power [mW]');
ylabel('HD output power [dBm]');
title('Noise Values from the Frequency Spectra');
set(gca,'Ytick',-90:1:0);
grid on;
hold on;
% Plot scaling-fit
plot(powerLO, fitY);
% Plot electronic noise
electronicNoise = ones(length(powerLO),1) * noise(1);
plot(powerLO,electronicNoise);
hold off;

% Saving figure
print(outputFilename,outputFiletype);

end

