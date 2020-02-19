function [] = g22ChvsDelayPhaseAveraged(channela, channelb)
% this function computes the second order correlation function between two channels 
% for a delay series if both channels dont possess a defined phase.

%% create data structs
dataStruct = struct('filename',{},'Delay', {},'NPhotonsCh1', {},'G2Ch1',{},'NPhotonsCh2',{},'G2Ch2',{},'G2Ch1Ch2',{});
dataStructLOonly = struct('filename',{},'number',{});

% %read delays from excel file
% delays = xlsread('2017-11-23-results','B2:B27');
% 
% for i = 1:size(delays)
%     dataStruct(i).Delay = delays(i);
% end

%% Create data overview
dispstat('','init','notquiet');
dispstat('Checking filenames ...','timestamp','keepthis','notquiet');
rawDataContents = dir('raw-data');

% LOwithDL-files
for name = {rawDataContents.name}
    % Loop only over *LOwithDL.raw files
    filename = cell2mat(name);
    if not(isempty(regexpi(filename,'.stamp','match'))) || not(isempty(regexpi(filename,'.cfg','match')))...
            ||(isempty(regexpi(filename,'LOwithDL','match')))
        continue
    end
    
    % Fetch number of measurement
    numberToken = regexpi(filename,'^([0123456789]*)-','tokens');
    number = str2double(cell2mat(numberToken{1}));
    dataStruct(number).filename = filename;
    
    % Get Delay
    %delayToken = regexpi(filename,'DL-([0123456789,]*)','tokens');
    delayToken = regexpi(filename,'raw([0123456789,]*)','tokens');
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
%for number = 1:2

    filenameSIG = dataStruct(number).filename;
    if isempty(filenameSIG)
        continue
    end
    
        % Use available quadrature dataset, if possible
    mainContents = dir();
    datasetExisting = 0;
    for name = {mainContents.name}
        mainName = cell2mat(name);
        if not(isempty(regexpi(mainName,['quadratureDataset-' strrep(num2str(filenameSIG),...
                '.raw','')],'match')))
            datasetExisting = 1;
            M = load(mainName);
            X1 = M.X1;
            X2 = M.X2;
            X3 = M.X3;
        else
            continue
        end
    end
    
    if datasetExisting == 0
        %find adequate LOonly-file
        LOnumber = max(LOnumbers(LOnumbers<=number));
        filenameLO = dataStructLOonly(LOnumber).filename;
         dispstat(['PrepareData number ' num2str(number)],...
            'timestamp','keepthis','notquiet');
        %prepare data
        [X1, X2, X3, ~,~] = prepareData(filenameLO, filenameSIG, 'Channels',[channela channelb],'Offset','local','Piezo','yes');        
%         [X1,X2,X3,~] = prepare3ChData(filenameLO, filenameSIG);
        save(['quadratureDataset-' strrep(num2str(filenameSIG),'.raw','.mat')], 'X1','X2','X3');
        %save quadratures
    end
    
    % First channel
     switch channela
        case 1
            Xa = X1;
        case 2
            Xa = X2;
        case 3
            Xa = X3;
    end
    Xa = Xa(:);
    [g2values,nPhotons] = g2(Xa,length(Xa));

    disp(['Delay: ',num2str(dataStruct(number).Delay)]) 
    disp(['First g2 of Channel 1: ',num2str(g2values)])
    disp(['First nPhotons: ',num2str(nPhotons)])
    dataStruct(number).G2Ch1 = g2values;
    dataStruct(number).NPhotonsCh1 = nPhotons;
    %save(strcat(num2str(i+1,'%02d'),'-LOwithSIG-DL',num2str(Delays(i)),'-Ch1.mat'),'X1','piezoSign1','theta1','uniformX1','uniformTheta1');

    % second channel
     switch channelb
        case 1
            Xb = X1;
        case 2
            Xb = X2;
        case 3
            Xb = X3;
    end
    Xb = Xb(:);
    [g2values,nPhotons] = g2(Xb,length(Xb));
    disp(['Delay: ',num2str(dataStruct(number).Delay)]) 
    disp(['First g2 of Channel 2: ',num2str(g2values)])
    disp(['First nPhotons: ',num2str(nPhotons)])
    dataStruct(number).G2Ch2 = g2values;
    dataStruct(number).NPhotonsCh2 = nPhotons;
    %save(strcat(num2str(i+1,'%02d'),'-LOwithSIG-DL-',num2str(Delays(i)),'-Ch2.mat'),'X2','piezoSign2','theta2','uniformX2','uniformTheta2');

    % correlation between both channels
    Xab = Xa(:)+Xb(:);
    g2values = g22Ch(Xa, Xb, Xab,length(Xa),length(Xb),length(Xab));
    disp(['Delay: ',num2str(dataStruct(number).Delay)])
    disp(['First g2 for two channels: ',num2str(g2values)])
    dataStruct(number).G2Ch1Ch2 = g2values;
   
end

%% Get Arrays with values

Delays = cell2mat({dataStruct.Delay});
NPhvsDLCh1 = cell2mat({dataStruct.NPhotonsCh1});
G2vsDLCh1 = cell2mat({dataStruct.G2Ch1});
NPhvsDLCh2 = cell2mat({dataStruct.NPhotonsCh2});
G2vsDLCh2 = cell2mat({dataStruct.G2Ch2});
G2vsDLCh1Ch2 = cell2mat({dataStruct.G2Ch1Ch2});

save(['G2-Ch1Ch2-channel1-' num2str(channela) '-channel2-' num2str(channelb) '.mat'],...
    'Delays','NPhvsDLCh1','G2vsDLCh1','NPhvsDLCh2','G2vsDLCh2','G2vsDLCh1Ch2');

% figure(1),plot((Delays-Delays(1))*10/3/1*2,NPhvsDLCh1,(Delays-Delays(1))*10/3/1*2,NPhvsDLCh2);
% figure(2),plot((Delays-Delays(1))*10/3/1*2,G2vsDLCh1,(Delays-Delays(1))*10/3/1*2,G2vsDLCh2);
% figure(3),plot((Delays-Delays(1))*10/3/1*2,G2vsDLCh1Ch2);

end
