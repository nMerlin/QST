function [] = g22ChDUvsDelay(period1, period2)
% % this function computes the second order correlation function between two channels 
% for a delay series if both channels possess a defined phase, i. e. are coherent.
% period1, period2: For each channel, the number of periods per piezo segment should be set
% according to the piezo modulation. For peak detection it is important to know how many wavelengths are
% located in one measured piezo segment. Optional: Implement automatic
% computation from config.
%E.g.: channel one 2 periods, channel two 4 periods.
% this function uses "deep uniforming" (uniformSampling2Channel.m) to
% ensure both channels have all phases uniformed AND their relative phase uniformed.  

%% create data structs
dataStruct = struct('filename',{},'Delay', {},'NPhotonsCh1', {},'G2Ch1',{},'NPhotonsCh2',{},'G2Ch2',{},'G2Ch1Ch2',{},'G2duCh1',{},'G2duCh2',{},'NPhduCh1',{},'NPhduCh2',{},'G2duCh1Ch2',{});
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

    filenameSIG = dataStruct(number).filename;
    if isempty(filenameSIG)
        continue
    end
    
    %find adequate LOonly-file
    LOnumber = max(LOnumbers(LOnumbers<=number));
    filenameLO = dataStructLOonly(LOnumber).filename;
    
    
     dispstat(['PrepareData number ' num2str(number)],...
        'timestamp','keepthis','notquiet');
    
    % prepare data for nonrandom phase
    [X1, X2, ~, piezoSign, ~] = prepareData(filenameLO,filenameSIG,'Channels',1:2,'Offset','global','Piezo','yes');
    
    offset1 = (mean(max(max(max(X1,1))))+mean(min(min(min(X1,1)))))/2;
    X1 = X1 - offset1;
    offset2 = (mean(max(max(max(X2,1))))+mean(min(min(min(X2,1)))))/2;
    X2 = X2 - offset2;
    % First channel      
    theta1 = computePhase(X1,ones(size(X1)),piezoSign,'Period',period1);
    % g2 estimation
    nG2 = 1;
    [g2values,nPhotons] = deal(zeros(nG2,1));
    for iG2=1:nG2
        [uniformX1,uniformTheta1] = seriesUniformSampling(X1(:),theta1(:),'NBins',200);
        [g2values(iG2),nPhotons(iG2)] = g2(uniformX1,length(uniformX1));
    end
    disp(['Delay: ',num2str(dataStruct(number).Delay)]) 
    disp(['First g2 of Channel 1: ',num2str(g2values(1))])
    disp(['First nPhotons: ',num2str(nPhotons(1))])
    G2Ch1 = g2values(1);
    dataStruct(number).G2Ch1 = g2values(1);
    NPhotonsCh1 = nPhotons(1);
    dataStruct(number).NPhotonsCh1 = nPhotons(1);
    %save(strcat(num2str(i+1,'%02d'),'-LOwithSIG-DL',num2str(Delays(i)),'-Ch1.mat'),'X1','piezoSign1','theta1','uniformX1','uniformTheta1');

    % second channel
    theta2 = computePhase(X2,ones(size(X2)),piezoSign,'Period',period2); 
    % g2 estimation
    nG2 = 1;
    [g2values,nPhotons] = deal(zeros(nG2,1));
    for iG2=1:nG2
        [uniformX2,uniformTheta2] = seriesUniformSampling(X2(:),theta2(:),'NBins',200);
        [g2values(iG2),nPhotons(iG2)] = g2(uniformX2,length(uniformX2));
    end
    disp(['Delay: ',num2str(dataStruct(number).Delay)]) 
    disp(['First g2 of Channel 2: ',num2str(g2values(1))])
    disp(['First nPhotons: ',num2str(nPhotons(1))])
    G2Ch2 = g2values(1);
    dataStruct(number).G2Ch2 = g2values(1);
    NPhotonsCh2 = nPhotons(1);
    dataStruct(number).NPhotonsCh2 = nPhotons(1);
    %save(strcat(num2str(i+1,'%02d'),'-LOwithSIG-DL-',num2str(Delays(i)),'-Ch2.mat'),'X2','piezoSign2','theta2','uniformX2','uniformTheta2');
    
    %correlation between both channels
    % g2 estimation
    theta12 = mod(theta2 - theta1,2*pi);
    nG2 = 1;
    [g2values] = deal(zeros(nG2,1));
    for iG2=1:nG2
        [uniformX12,uniformTheta12] = seriesUniformSampling(X1(:)+X2(:),theta12(:),'NBins',200);
        g2values(iG2) = g22Ch(uniformX1, uniformX2, uniformX12,length(uniformX1),length(uniformX2),length(uniformX12));
    end
    disp(['Delay: ',num2str(dataStruct(number).Delay)])
    disp(['First g2 for two channels: ',num2str(g2values(1))])
    G2Ch1Ch2 = g2values(1);
    dataStruct(number).G2Ch1Ch2 = g2values(1);
    
    %deep uniform
    [uniform2ChX1,uniform2ChX2,uniform2ChTheta1,uniform2ChTheta2] = uniformSampling2Channel(X1,X2,theta1,theta2,2);
    g2values12sum = g22Ch(uniform2ChX1, uniform2ChX2, uniform2ChX1+uniform2ChX2,length(uniform2ChX1),length(uniform2ChX2),length(uniform2ChX1+uniform2ChX2));
    [g2values2Ch1,nPhotons2Ch1] = g2(uniform2ChX1,length(uniform2ChX1));
    [g2values2Ch2,nPhotons2Ch2] = g2(uniform2ChX2,length(uniform2ChX2));
    %disp(['Delay: ',num2str(dataStruct(number).Delay)]) 
    disp(['First g2 for two channel deep uniform: ',num2str(g2values12sum(1))])
    %disp(['First nPhotons: ',num2str(nPhotons(1))])
    %G2duCh1 = g2values2Ch1(1);
    dataStruct(number).G2duCh1 = g2values2Ch1(1);
    dataStruct(number).G2duCh2 = g2values2Ch2(1);
    dataStruct(number).NPhduCh1 = nPhotons2Ch1(1);
    dataStruct(number).NPhduCh2 = nPhotons2Ch2(1);
    dataStruct(number).G2duCh1Ch2 = g2values12sum(1);
    
    
    save(strcat(filenameSIG, '.mat'),'X1','X2','piezoSign','theta1','theta2','uniformX1','uniformX2','uniformTheta1','uniformTheta2',...
    'theta12','G2Ch1','G2Ch2','NPhotonsCh1','NPhotonsCh2','G2Ch1Ch2','g2values12sum','g2values2Ch1','g2values2Ch2','nPhotons2Ch1','nPhotons2Ch2');
   
end


%% Get Arrays with values

Delays = cell2mat({dataStruct.Delay});
NPhvsDLCh1 = cell2mat({dataStruct.NPhotonsCh1});
G2vsDLCh1 = cell2mat({dataStruct.G2Ch1});
NPhvsDLCh2 = cell2mat({dataStruct.NPhotonsCh2});
G2vsDLCh2 = cell2mat({dataStruct.G2Ch2});
G2vsDLCh1Ch2 = cell2mat({dataStruct.G2Ch1Ch2});
G2vsDLduCh1 = cell2mat({dataStruct.G2duCh1});
G2vsDLduCh2 = cell2mat({dataStruct.G2duCh2});
NPhvsDLduCh1 = cell2mat({dataStruct.NPhduCh1});
NPhvsDLduCh2 = cell2mat({dataStruct.NPhduCh2});
G2vsDLduCh1Ch2 = cell2mat({dataStruct.G2duCh1Ch2});


save(['G2-Ch1Ch2-period1-' num2str(period1) '-period2-' num2str(period2) '.mat'],...
    'Delays','NPhvsDLCh1','G2vsDLCh1','NPhvsDLCh2','G2vsDLCh2','G2vsDLCh1Ch2','G2vsDLduCh1','G2vsDLduCh2','NPhvsDLduCh1','NPhvsDLduCh2','G2vsDLduCh1Ch2');

figure(1),plot((Delays-Delays(1))*10/3/1*2,NPhvsDLCh1,(Delays-Delays(1))*10/3/1*2,NPhvsDLCh2);
figure(2),plot((Delays-Delays(1))*10/3/1*2,G2vsDLCh1,(Delays-Delays(1))*10/3/1*2,G2vsDLCh2);
figure(3),plot((Delays-Delays(1))*10/3/1*2,G2vsDLduCh1Ch2);

end
