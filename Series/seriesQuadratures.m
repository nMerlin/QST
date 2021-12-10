function [] = seriesQuadratures(varargin)
% % this function computes the quadratures for one channel 
% for a series.
% 'Channels':  Array with the numbers of the channels that should be
%   evaluated.
%E.g.: 2 periods for 50 Hz modulation and 2 um amplitude 
%   'CorrRemove': Array that should be as long as the 'Channels' array, to decide for each channel.
%   When true, the correlations of quadratures with precedent
%   quadratures are removed. This may take a few minutes. Recommended only for random phase. 
%   'Offset': Array that should be as long as the 'Channels' array, to decide for each channel. 
%   when 'global', a global offset is subtracted from the
%   quadratures. Recommended when phase is not random. When 'local', offset
%   is removed piecewise. Recommended for random phase. 

%% Validate and parse input arguments
p = inputParser;
defaultChannels = 1;
addParameter(p,'Channels',defaultChannels,@isnumeric);
defaultOffset = ["local" "local" "local"]; % vs "global"
addParameter(p,'Offset',defaultOffset,@isstring);
defaultPiezo = 'yes'; 
addParameter(p,'Piezo',defaultPiezo);
defaultCorrRemove = [true,true,true];
addParameter(p,'CorrRemove',defaultCorrRemove,@islogical);
parse(p,varargin{:});
c = struct2cell(p.Results);
[channels,corrRemove,offset,piezo] = c{:};

%% create data structs
%dataStruct = struct('filename',{},'Delay', {},'nPhMean', {},'G2',{});
dataStruct = struct('filename',{});
dataStructLOonly = struct('filename',{},'number',{});

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
   
end

% LOonly-files
for name = {rawDataContents.name}
    % Loop only over *LOonly.raw files
    filename = cell2mat(name);
    if not(isempty(regexpi(filename,'.stamp','match'))) || not(isempty(regexpi(filename,'.cfg','match')))...
            || isempty(regexpi(filename,'LOonly','match'))
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
    %LOnumber = min(LOnumbers);
    filenameLO = dataStructLOonly(LOnumber).filename;
    
    
     dispstat(['PrepareData number ' num2str(number)],...
        'timestamp','keepthis','notquiet');
    
    if ~exist('mat-data','dir')
    mkdir('mat-data')
    end
    
    offsetname = char(offset);
    offsetname = offsetname(:)';
    
    if strcmp(piezo,'yes')  
        [X1, X2, X3, piezoSign, ~] = prepareData(filenameLO,filenameSIG,...
            'Channels',channels,'Offset',offset,'Piezo',piezo,'CorrRemove',corrRemove);        
        cd('mat-data');
        save(strcat(filenameSIG,'-offset-',offsetname,'-piezo-',piezo,'-corrRemove-',num2str(corrRemove),'.mat'),'X1','X2','X3','piezoSign');
        cd('..');

    else
       [X1, X2, X3, ~, ~] = prepareData(filenameLO,filenameSIG,'Channels',...
           channels,'Offset',offset,'Piezo',piezo,'CorrRemove',corrRemove);   
        cd('mat-data');
        save(strcat(filenameSIG,'-offset-',offsetname,'-piezo-',piezo,'-corrRemove-',num2str(corrRemove), '.mat'),'X1','X2','X3');
        cd('..');
    end
     
end

end
