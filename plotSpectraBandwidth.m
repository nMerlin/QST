function [  ] = plotSpectraBandwidth( varargin )
%PLOTPOWERSPECTRA Summary of this function goes here
%   Detailed explanation goes here
%   Assumes the datafiles in the directory 'raw-data'. The filenames should
%   be in the format '01-10mW-FFT.txt'.
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
dataStruct = struct('filename',{},'powerLO',{},'BW',{});
fcount = 0;
outputFilename1 = 'Compare-power-spectra_smoothed.jpg';
outputFiletype1 = '-djpeg';
outputFilename2 = 'Bandwith.jpg';
outputFiletype2 = '-djpeg';
AveragingLimits = [5 40]; % in MHz
span= 5; %Span of the moving average for smoothing of spectra
 
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
    if isempty(regexpi(filename,'FFT.txt','match'))
        continue
    end
    
    fcount = fcount + 1;
    dataStruct(fcount).filename = filename;
    
    % Get LO power
    powerToken = regexpi(filename,'-([0123456789.]*)mW','tokens');
    dataStruct(fcount).powerLO = str2double(cell2mat(powerToken{1}));
   
end

dispstat('Loading data for plots ...','timestamp','keepthis',quiet);
%%% Process data for plotting
powerLO = cell2mat({dataStruct.powerLO});

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

% Plot the spectra
dispstat('Plotting', 'timestamp','keepthis',quiet);
close all;

startIndex = find(frequencyAxis > AveragingLimits(1)*1000000,1);
stopIndex = find(frequencyAxis > AveragingLimits(2)*1000000,1);

for k=1:size(dataStruct,2)
    fft = dlmread(dataStruct(k).filename,' ',12,0);
    fftsmooth = smooth(fft(:,2), span);
    fft(:,1) = frequencyAxis;
   
    Average = mean(fftsmooth(startIndex:stopIndex));
    IndexBW = find(fftsmooth<(Average - 3),1);
    dataStruct(k).BW = frequencyAxis(IndexBW);
    
    plot(fft(:,1)/1000000, fftsmooth,'DisplayName', strcat(num2str(dataStruct(k).powerLO), ' mW '));
    hold on;
end
cd('..');
    
axis([xMin xMax yMin yMax]);
xlabel('Frequency [MHz]');
ylabel('HD output power [dBm]');
legend('Location','northeast');
hold off;
% Saving figure
print(outputFilename1,outputFiletype1);

BW = cell2mat({dataStruct.BW});

plot(powerLO,BW/1000000,'o');
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultTextInterpreter','latex');
xlabel('P_{LO} [mW]');
ylabel('Bandwidth [MHz]');
% Saving figure
print(outputFilename2,outputFiletype2);

dispstat('Finished!','timestamp',quiet);

end

