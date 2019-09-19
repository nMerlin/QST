function [] = seriesQuadratures(varargin)
% % this function computes the quadratures for one channel 
% for a series.
% 'Channels':  Array with the numbers of the channels that should be
%   evaluated.
% 'Period': The number of periods per piezo segment should be set
% according to the piezo modulation. For peak detection it is important to know how many wavelengths are
% located in one measured piezo segment. Optional: Implement automatic
% computation from config.
%E.g.: 2 periods for 50 Hz modulation and 2 um amplitude
% 'Phase': decide whether the phase is random or fixed. 
% 

%% Validate and parse input arguments
p = inputParser;
defaultChannels = 1;
addParameter(p,'Channels',defaultChannels,@isnumeric);
defaultPeriod1 = 2;
addParameter(p,'Period',defaultPeriod1,@isnumeric);
defaultPhase = 'random';
addParameter(p,'Phase',defaultPhase);
parse(p,varargin{:});
c = struct2cell(p.Results);
[channels,period1,phase] = c{:};

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
    if not(isempty(regexpi(filename,'.raw.','match')))...
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
    
    if strcmp(phase,'nonrandom')   
%        One channel with non random phase 
        [X1, X2, X3, piezoSign1, ~] = prepareData(filenameLO,filenameSIG,'Channels',channels,'Offset','global','Piezo','yes');
        if any(any(X1))
            X = X1;
        elseif any(any(X2))
            X = X2;
        elseif any(any(X3))
            X = X3;
        end
        
        theta1 = computePhase(X,ones(size(X)),piezoSign1,'Period',period1); 
        cd('Quadratures');
        save(strcat(filenameSIG, '.mat'),'X','piezoSign1','theta1');
        cd('..');

    elseif strcmp(phase,'random')
   % One channel with random phase 
       [X1, X2, X3, ~, ~] = prepareData(filenameLO,filenameSIG,'Channels',channels,'Offset','local','Piezo','no');
            if any(any(X1))
                X = X1;
            elseif any(any(X2))
                X = X2;
            elseif any(any(X3))
                X = X3;
            end
   
       
         cd('Quadratures');
        save(strcat(filenameSIG, '.mat'),'X');
        cd('..');
    end
     
end

end
