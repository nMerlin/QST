function [] = g21ChvsDelayFromMat(varargin)
% % this function computes the second order correlation function between two channels 
% for a delay series if both channels possess a defined phase, i. e. are coherent.
% period1, period2: For each channel, the number of periods per piezo segment should be set
% according to the piezo modulation. For peak detection it is important to know how many wavelengths are
% located in one measured piezo segment. Optional: Implement automatic
% computation from config.
%E.g.: mostly channel one 2 periods

%% Validate and parse input arguments

%% create data structs
dataStruct = struct('filename',{},'Delay', {},'NPhotons', {},'G2',{});

%% Create data overview
dispstat('','init','notquiet');
dispstat('Checking filenames ...','timestamp','keepthis','notquiet');
rawDataContents = dir();

% LOwithDL-files
for name = {rawDataContents.name}
    % Loop only over *LOwithDL.raw files
    filename = cell2mat(name);
    if (isempty(regexpi(filename,'LOwithDL.raw.mat','match')))
        continue
    end
    
    % Fetch number of measurement
    numberToken = regexpi(filename,'^([0123456789]*)-','tokens');
    number = str2double(cell2mat(numberToken{1}));
    dataStruct(number).filename = filename;
    
    % Get Delay
    %delayToken = regexpi(filename,'DL-([0123456789,]*)','tokens');
    delayToken = regexpi(filename,'-([0123456789,]*)mm','tokens');
    dataStruct(number).Delay = str2double(strrep(cell2mat(delayToken{1}),',','.'));
    %dataStruct(number).Delay = str2double(cell2mat(delayToken{1}));
end

%% Loop over files

for number = 1:size(dataStruct,2)
%for number = 1:2

    filenameSIG = dataStruct(number).filename;
    if isempty(filenameSIG)
        continue
    end
    T = load(filenameSIG); 
 
    dataStruct(number).NPhotons = T.nPhotons;
    dataStruct(number).G2 = T.g2values;
    
end

%% Get Arrays with values

Delays = cell2mat({dataStruct.Delay});
NPhvsDLCh1 = cell2mat({dataStruct.NPhotons});
G2vsDLCh1 = cell2mat({dataStruct.G2});

save('G2-1Ch.mat','Delays','NPhvsDLCh1','G2vsDLCh1');

figure(1),plot((Delays-Delays(1))*10/3/1*2,NPhvsDLCh1);
figure(2),plot((Delays-Delays(1))*10/3/1*2,G2vsDLCh1);


end
