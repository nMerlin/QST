function [X1, X2, X3, piezoSign, configSIG] = prepareData(filenameLO, filenameSIG, varargin)
%PREPARE3CHDATA Returns quadratures of a 3-Channel-Measurement
%
%   filenameLO - filename of the LO-data used for correct normalization
%   filenameSIG - filename of the raw 3-Channel-Data
%
%   Options:
%   'Channels': Array with the numbers of the channels that should be
%   evaluated.
%   'CorrRemove': Array that should be as long as the 'Channels' array, to decide for each channel.
%   When true, the correlations of quadratures with precedent
%   quadratures are removed. This may take a few minutes. Recommended only for random phase. 
%   'Offset': Array that should be as long as the 'Channels' array, to decide for each channel. 
%   when 'global', a global offset is subtracted from the
%   quadratures. Recommended when phase is not random. When 'local', offset
%   is removed piecewise. Recommended for random phase. 
%   'Piezo': when 'yes', the quadrature matrices are cut into piezo segments.
%   'Save': when 'yes', the quadrature data is saved.
%
% Output Arguments:
%   piezoSign: +1 means piezo moves in the first segment from 0 to 2 um
%              -1 means piezo moves in the first segment from 2 to 0 um
%               0 if no piezo modulation used

%% Validate and parse input arguments
p = inputParser;
defaultChannels = 1:3;
addParameter(p,'Channels',defaultChannels,@isnumeric);
defaultPickingFactor = 1;
addParameter(p,'PickingFactor',defaultPickingFactor,@isnumeric);
defaultCorrRemove = [true,true,true];
addParameter(p,'CorrRemove',defaultCorrRemove,@islogical);
defaultDutyCycle = 1/3; % Integration Duty Cycle
addParameter(p,'DutyCycle',defaultDutyCycle,@isnumeric);
defaultSave = 'no'; % save quadratures as mat
addParameter(p,'Save',defaultSave);
defaultOffset = ["local" "local" "local"]; % vs 'global'
addParameter(p,'Offset',defaultOffset,@isstring);
defaultPiezo = 'yes'; 
addParameter(p,'Piezo',defaultPiezo);
parse(p,varargin{:});
c = struct2cell(p.Results);
[channels,corrRemove,dutycycle,offset,pickingFactor,piezo,saveOption] = c{:};


%% preliminary values
[X1,X2,X3, piezoSign] = deal(zeros(1));
maxChannelNumber = 3; % maximum number of channels

%% load data
CALIBRATION_CH1 = 4.596047840078126e-05; % Ampere per Volt

%Norm ist the factor in the relation between the quadratures and the ladder
%operators: q = Norm*(a^{+} + a), p = Norm*i*(a^{+} - a)
%typical values are 1/sqrt(2) or 1/2.
Norm = 1/sqrt(2);

dispstat('','init','timestamp','keepthis',0);
dispstat('Load LO data','timestamp','keepthis',0);
[data8bitLO,configLO,~]=load8BitBinary(filenameLO,'dontsave');

dispstat('Load Signal data','timestamp','keepthis',0);
[data8bitSIG,configSIG,timestamps]=load8BitBinary(filenameSIG,'dontsave');

%% Compute number of LO photons
dispstat('Computing number of LO photons for Channels.', ...
    'timestamp','keepthis',0);
if pickingFactor == 20
    XLO = computeQuadratures(data8bitLO(:,:,channels),configLO, ...
        CALIBRATION_CH1,'LocationOffset',8, ...
        'IntegrationWindow',30,'MinPeakDistance',75,'DutyCycle',dutycycle);
else
    XLO = computeQuadratures(data8bitLO(:,:,channels),configLO, ...
    CALIBRATION_CH1,'DutyCycle',dutycycle);
end

% fill the uncomputed channels with zeros.
Xfull = zeros(size(XLO,1),size(XLO,2),maxChannelNumber);
Xfull(:,:,channels) = XLO;
XLO = Xfull;

if any(corrRemove)
        %Removing correlations with precedent pulses for the LO
        XLO(:,:,channels) = bsxfun(@minus, XLO(:,:,channels), mean(XLO(:,:,channels)));
        dispstat('Removing correlations from LO... ','timestamp','keepthis',0);
        XLO = correlationCompensation(XLO);
        XLO = XLO(1:end-1,:,:); % remove linearly dependent last vector
end

% Calculate the variance piece-wise to compensate slow drifts (e.g.
% piezos)
NLO = mean(var(XLO));

%% Compute quadratures for target quantum state
dispstat('Computing target quadratures for Channels', ...
    'timestamp','keepthis',0);
if pickingFactor == 20
    X = computeQuadratures(data8bitSIG(:,:,channels),configSIG, ...
        CALIBRATION_CH1,'LocationOffset',8, ...
        'IntegrationWindow',30,'MinPeakDistance',75,'DutyCycle',dutycycle);
else
    X = computeQuadratures(data8bitSIG(:,:,channels),configSIG, ...
        CALIBRATION_CH1,'DutyCycle',dutycycle);
end

% fill the uncomputed channels with zeros.
Xfull = zeros(size(X,1),size(X,2),maxChannelNumber);
Xfull(:,:,channels) = X;
X = Xfull;

%% calibration to vacuum state and offset removal
for iCh = channels
    % Calibration of quadratures to vacuum state
    X(:,:,iCh) = Norm * X(:,:,iCh) / sqrt(NLO(iCh));
    
    % Removing Offsets
    if strcmp(offset(iCh),'local')  % remove piecewise offset
        dispstat(['Removing piecewise offset from channel ',num2str(iCh), ...
            ' ...'],'timestamp','keepthis',0);
        X(:,:,iCh) = bsxfun(@minus, X(:,:,iCh), mean(X(:,:,iCh)));
        Xrem = X(:,:,iCh);
        
        if corrRemove(iCh) %Removing correlations with precedent pulses
            dispstat(['Removing correlations from channel ',num2str(iCh), ...
                ' ...'],'timestamp','keepthis',0);
            X(:,:,iCh) = correlationCompensation(X(:,:,iCh));
            Xrem = X(1:end-1,:,iCh);
        end 
             
    elseif strcmp(offset(iCh),'global') % remove global offset
        dispstat(['Removing global offset from channel ',num2str(iCh), ...
            ' ...'],'timestamp','keepthis',0);
        Xrem = X(:,:,iCh) - mean(mean(X(:,:,iCh)));       
    end
    
%% Cut the raw data into segments of equal length according to piezo, if piezo modulation was used
    % modulation
    if strcmp(piezo,'yes')
        dispstat(['Reshaping channel ',num2str(iCh), ...
            ' into piezo segments ...'],'timestamp','keepthis',0);
        [Xrem,piezoSign] = piezoSegments(timestamps,Xrem,'cut');       
        if strcmp(offset(iCh),'global')
            %compute offset from mean maxima and minima to prevent they are
            %biased
            biase = (mean(max(max(max(Xrem,1))))+mean(min(min(min(Xrem,1)))))/2;
            Xrem = Xrem - biase;
        end
    end
     
%% asign channels
    switch iCh
        case 1
            X1 = Xrem;
        case 2
            X2 = Xrem;
        case 3
            X3 = Xrem;
    end
  
end % iCh

%% Saving important Variables (to delete the raw data manually)
if strcmp(saveOption,'yes')
    dispstat('Saving ...','timestamp','keepthis',0);
    dateString = datestr(datetime('now'),'yyyymmddTHHMMSS');
    save(['quadratureDataset-', strrep(num2str(filenameSIG),...
        '.raw','') dateString '.mat'], 'X1','X2','X3','piezoSign');
end


end % prepareData

