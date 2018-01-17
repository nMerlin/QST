function [  ] = plotPowerSpectra( varargin )
%PLOTPOWERSPECTRA Summary of this function goes here
%   Detailed explanation goes here
%   Assumes the datafiles in the directory 'raw-data'. The filenames should
%   be in the format '01-10mW-filtered-FFT.txt' or '01-10mW-pure-FFT.txt'.
%   The files need to contain a header created by SBench.

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
dataStruct = struct('filename',{},'powerLO',{},'filtered',{},'singlediode',{});
fcount = 0;
outputFilename = 'power-spectra.jpg';
outputFiletype = '-djpeg';
clearanceAveragingLimits = [20 40]; % in MHz
xMin = 0;
xMax = 200;
yMin = -90;
yMax = 0;

%%% Create data overview
dispstat('','init',quiet);
dispstat('Checking filenames ...','timestamp','keepthis',quiet);
rawDataContents = dir('raw-data');
for name = {rawDataContents.name}
    % Loop only over *-FFT.txt files & get filename
    filename = cell2mat(name);
    if isempty(regexpi(filename,'FFT','match'))
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

% Load the three datasets of interest
dispstat('Load the three datasets of interest ...', ...
    'timestamp','keepthis',quiet);
for k=1:size(dataStruct,2)
    if dataStruct(k).powerLO == maxPower && dataStruct(k).filtered == 0
        response = dlmread(dataStruct(k).filename,' ',12,0);
        response(:,1) = frequencyAxis;
    elseif dataStruct(k).powerLO == maxPower && dataStruct(k).filtered == 1
        filteredResponse = dlmread(dataStruct(k).filename,' ',12,0);
        filteredResponse(:,1) = frequencyAxis;
    elseif dataStruct(k).powerLO == minPower && dataStruct(k).filtered == 0
        electronicNoise = dlmread(dataStruct(k).filename,' ',12,0);
        electronicNoise(:,1) = frequencyAxis;
    elseif dataStruct(k).powerLO == 0.1 && dataStruct(k).filtered == 0 ...
            && dataStruct(k).singlediode == 0
        cmrrTwoDiodes = dlmread(dataStruct(k).filename,' ',12,0);
        cmrrTwoDiodes(:,1) = frequencyAxis;
    elseif dataStruct(k).powerLO == 0.1 && dataStruct(k).filtered == 0 ...
            && dataStruct(k).singlediode == 1
        cmrrOneDiode = dlmread(dataStruct(k).filename,' ',12,0);
        cmrrOneDiode(:,1) = frequencyAxis;
    end
end
cd('..');

% Calculate the shot noise clearance
dispstat('Calculate shot noise clearance ...','timestamp','keepthis',quiet);
startIndex = find(frequencyAxis > clearanceAveragingLimits(1)*1000000,1);
stopIndex = find(frequencyAxis > clearanceAveragingLimits(2)*1000000,1);
responseAverage = mean(response(startIndex:stopIndex,2));
electronicNoiseAverage = mean(electronicNoise(startIndex:stopIndex,2));
clearance = responseAverage - electronicNoiseAverage;

% Calculate the CMRR
cmrrStart = find(frequencyAxis > 75*1000000,1);
cmrrStop = find(frequencyAxis > 76*1000000,1);
[peakTwoDiodes,~]=max(findpeaks(cmrrTwoDiodes(cmrrStart:cmrrStop,2)));
[peakOneDiode,~]=max(findpeaks(cmrrOneDiode(cmrrStart:cmrrStop,2)));

dispstat('Plotting ...','timestamp','keepthis',quiet);
%%% Plotting
close all;
% Plot data curves
plot(response(:,1)/1000000,response(:,2));
axis([xMin xMax yMin yMax]);
xlabel('Frequency [MHz]');
ylabel('HD output power [dBm]');
hold on;
plot(filteredResponse(:,1)/1000000,filteredResponse(:,2));
plot(electronicNoise(:,1)/1000000,electronicNoise(:,2));
hold off;
legend(strcat('Response (',num2str(maxPower,4),' mW)'), ...
    strcat('Filtered Response (',num2str(maxPower,4),' mW)'), ...
    strcat('Response (',num2str(minPower),' mW)'), ...
    'Location','northeast');

% Plot shot noise clearance
middleIndex = round((startIndex + stopIndex)/2);
xMid = frequencyAxis(middleIndex)/1000000;
h = annotation('doublearrow');
set(h,'parent',gca, ...
    'Position',[xMid electronicNoiseAverage 0 clearance]);
h = annotation('textbox');
set(h,'parent',gca, ...
    'Position',[xMid+2 electronicNoiseAverage+2*clearance/3 0 0], ...
    'String',strcat(num2str(clearance,3),' dBm'), ...
    'FitBoxToText','on','EdgeColor','none');

% Plot inset for CMRR
insetAx = axes('Parent',gcf,'Position',[0.22 0.6 0.15 0.25]);
plot(cmrrTwoDiodes(:,1)/1000000,cmrrTwoDiodes(:,2),'-k');
hold on;
plot(cmrrOneDiode(:,1)/1000000,cmrrOneDiode(:,2),'Color',[0.5 0.5 0.5]);
hold off;
set(insetAx,'XLim',[65 85], ...
    'YLim',[-90 -10], ...
    'FontSize',8);
xlabel('Frequency [MHz]');
ylabel('HD output power [dBm]');
title(['CMRR = ' num2str(peakOneDiode-peakTwoDiodes,3) ' dBm']);

% Saving figure
print(outputFilename,outputFiletype);

dispstat('Finished!','timestamp',quiet);

end

