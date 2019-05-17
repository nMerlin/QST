function [] = g21ChvsDelayFromMat(varargin)
% % this function computes the second order correlation function between two channels 
% for a delay series if both channels possess a defined phase, i. e. are coherent.
% period1, period2: For each channel, the number of periods per piezo segment should be set
% according to the piezo modulation. For peak detection it is important to know how many wavelengths are
% located in one measured piezo segment. Optional: Implement automatic
% computation from config.
%E.g.: mostly channel one 2 periods

%% Validate and parse input arguments
p = inputParser;
defaultPeriod1 = 2;
addParameter(p,'Period',defaultPeriod1,@isnumeric);
defaultPhase = 'random';
addParameter(p,'Phase',defaultPhase);
parse(p,varargin{:});
c = struct2cell(p.Results);
[period1,phase] = c{:};

%% create data structs
dataStruct = struct('filename',{},'Delay', {},'NPhotonsCh1', {},'G2Ch1',{});

%% Create data overview
dispstat('','init','notquiet');
dispstat('Checking filenames ...','timestamp','keepthis','notquiet');
rawDataContents = dir('mat-data')

% LOwithDL-files
for name = {rawDataContents.name}
    % Loop only over *LOwithDL.raw files
    filename = cell2mat(name);
    if not(isempty(regexpi(filename,'.mat','match')))...
            ||(isempty(regexpi(filename,'LOwithDL','match')))
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
    load(filenameSIG); 
    
    if strcmp(phase,'nonrandom')   
       % One channel with non random phase (change for random phase!)
        % g2 estimation
        nG2 = 1;
        [g2values,nPhotons] = deal(zeros(nG2,1));
        for iG2=1:nG
            [g2values(iG2),nPhotons(iG2)] = g2(uniformX1,length(uniformX1));
        end
        g2values = mean(g2values);
        nPhotons = mean(nPhotons);

    elseif strcmp(phase,'random')
   % One channel with random phase 
        % g2 estimation
        [g2values,nPhotons] = g2(X1,length(X1));
    end
    
    disp(['Delay: ',num2str(dataStruct(number).Delay)])   % test assess theta here!
    disp(['First g2 of Channel 1: ',num2str(g2values)])
    disp(['First nPhotons: ',num2str(nPhotons)])
    dataStruct(number).G2Ch1 = g2values; % test assess theta here!
    dataStruct(number).NPhotonsCh1 = nPhotons;
     
end

%% Get Arrays with values

Delays = cell2mat({dataStruct.Delay});
NPhvsDLCh1 = cell2mat({dataStruct.NPhotonsCh1});
G2vsDLCh1 = cell2mat({dataStruct.G2Ch1});

save(['G2-Ch1-period1-' num2str(period1) '.mat'],...
    'Delays','NPhvsDLCh1','G2vsDLCh1');
% figure(1),plot(NPhvsDLCh1);
% figure(2),plot(G2vsDLCh1);

%figure(1),plot((Delays-Delays(1))*10/3/1*2,NPhvsDLCh1);
%figure(2),plot((Delays-Delays(1))*10/3/1*2,G2vsDLCh1);

end
