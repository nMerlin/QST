function [] = g21ChvsDelay(varargin)
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
dataStruct = struct('filename',{},'Delay', {},'NPhotons', {},'G2',{});
dataStructLOonly = struct('filename',{},'number',{});

%% Create data overview
dispstat('','init','notquiet');
dispstat('Checking filenames ...','timestamp','keepthis','notquiet');
rawDataContents = dir('raw-data');

% LOwithDL-files
for name = {rawDataContents.name}
    % Loop only over *LOwithDL.raw files
    filename = cell2mat(name);
    if not(isempty(regexpi(filename,'.raw.','match')))...
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

% LOonly-files
for name = {rawDataContents.name}
    % Loop only over *LOonly.raw files
    filename = cell2mat(name);
    if not(isempty(regexpi(filename,'LOonly.raw.','match')))...
            || isempty(regexpi(filename,'LOonly.raw','match'))
        continue
    end
    
    % Fetch number of measurement
    numberToken = regexpi(filename,'^([0123456789]*)-','tokens');
    number = str2double(cell2mat(numberToken{1}));
    dataStructLOonly(number).filename = filename;
    dataStructLOonly(number).number = number;
end
LOnumbers = cell2mat({dataStructLOonly.number});

%% Loop over files

for number = 1:size(dataStruct,2)
%for number = 2:34
    
    filenameSIG = dataStruct(number).filename;
    if isempty(filenameSIG)
        continue
    end
    
    %find adequate LOonly-file
    LOnumber = max(LOnumbers(LOnumbers<=number));
    filenameLO = dataStructLOonly(LOnumber).filename;
    
    
     dispstat(['PrepareData number ' num2str(number)],...
        'timestamp','keepthis','notquiet');
    
    if strcmp(phase,'nonrandom')   
%        One channel with non random phase (change for random phase!)
        [~, X2, ~, piezoSign1, ~] = prepareData(filenameLO,filenameSIG,'Channels',2,'Offset','global','Piezo','yes');
        X1 = X2;
        theta1 = computePhase(X1,ones(size(X1)),piezoSign1,'Period',period1); 
%         %g2 estimation
%         load([filenameSIG '.mat']);
        nG2 = 1;
        [g2values,nPhotons] = deal(zeros(nG2,1));
        for iG2=1:nG2
            [uniformX1,uniformTheta1] = seriesUniformSampling(X1(:),theta1(:),'NBins',100);
            [g2values(iG2),nPhotons(iG2)] = g2(uniformX1,length(uniformX1));
        end
        g2values = mean(g2values);
        nPhotons = mean(nPhotons);
        save(strcat(filenameSIG, '.mat'),'X1','piezoSign1','theta1','uniformX1','uniformTheta1','g2values','nPhotons');

    elseif strcmp(phase,'random')
   % One channel with random phase 
        [~, X2, ~, piezoSign1, ~] = prepareData(filenameLO,filenameSIG,'Channels',2,'Offset','local','Piezo','yes');
        X1 = X2;
        X1 = X1(:);
        % g2 estimation
        [g2values,nPhotons] = g2(X1,length(X1));
        save(strcat(filenameSIG, '.mat'),'X1','piezoSign1','g2values','nPhotons');
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

figure(1),plot((Delays-Delays(1))*10/3/1*2,NPhvsDLCh1);
figure(2),plot((Delays-Delays(1))*10/3/1*2,G2vsDLCh1);

end
