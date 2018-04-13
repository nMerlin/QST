function [  ] = plotAnalog( varargin )
%PLOTPOWERSPECTRA Summary of this function goes here
%   Detailed explanation goes here
%   Assumes the datafiles in the directory 'raw-data'. The filenames should
%   be in the format '01-10mW-Analog.txt'.
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
dataStruct = struct('filename',{},'powerLO',{});
fcount = 0;
outputFilename1 = 'Analog';
outputFiletype1 = '-djpeg';

xMin = 0;
xMax = 60;
yMin = -1;
yMax = 1;

%%% Create data overview
dispstat('','init',quiet);
dispstat('Checking filenames ...','timestamp','keepthis',quiet);
rawDataContents = dir('raw-data');
for name = {rawDataContents.name}
    % Loop only over *-Analog.txt files & get filename
    filename = cell2mat(name);
    if isempty(regexpi(filename,'Analog.txt','match'))
        continue
    end
    
    fcount = fcount + 1;
    dataStruct(fcount).filename = filename;
    
    % Get LO power
    powerToken = regexpi(filename,'-([0123456789.]*)mW','tokens');
    dataStruct(fcount).powerLO = str2double(cell2mat(powerToken{1}));
   
end

cd('raw-data');

% Plot 
dispstat('Plotting', 'timestamp','keepthis',quiet);
close all;

for k=1:size(dataStruct,2)
    Analog = dlmread(dataStruct(k).filename,' ',12,0);
    Analog(:,1) = (Analog(:,1)-Analog(1,1))*10^9; 
    
    plot(Analog(:,1), Analog(:,2),'DisplayName', strcat(num2str(dataStruct(k).powerLO), ' mW '));
    hold on;
end
cd('..');
    
axis([xMin xMax yMin yMax]);
%xlabel('time [ns]');
%ylabel('HD output voltage [V]');
xlabel('Zeit (ns)');
ylabel('Signal (V)');
legend('Location','northeast');
hold off;
% Saving figure
%print([outputFilename1,'.jpg'],outputFiletype1);
saveA5Landscape(outputFilename1);

dispstat('Finished!','timestamp',quiet);

end

