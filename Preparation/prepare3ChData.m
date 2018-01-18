function [X1, X2, X3, piezoSign] = prepare3ChData(filenameLO, filenameSIG, varargin)
%PREPARE3CHDATA Returns quadratures of a 3-Channel-Measurement
%
% The quadrature matrices are already cut into piezo segments.
%
%   filenameLO - filename of the LO-data used for correct normalization
%   filenameSIG - filename of the raw 3-Channel-Data
%
% With 'CorrRemove','yes', you can have the correlations of quadratures
% with precedent quadratures removed. This may take a few minutes.
%
% Output Arguments:
%   piezoSign: +1 means piezo moves in the first segment from 0 to 2 um
%              -1 means piezo moves in the first segment from 2 to 0 um

%% Validate and parse input arguments
p = inputParser;
defaultPickingFactor = 1;
addParameter(p,'PickingFactor',defaultPickingFactor,@isnumeric);
defaultCorrRemove = 'no';
addParameter(p,'CorrRemove',defaultCorrRemove);
defaultDutyCycle = 1/3; % Integration Duty Cycle
addParameter(p,'DutyCycle',defaultDutyCycle,@isnumeric);
parse(p,varargin{:});
c = struct2cell(p.Results);
[corrRemove,dutycycle,pickingFactor] = c{:};

%%

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

% Compute number of LO photons
dispstat('Computing number of LO photons for Channels 1 to 3.', ...
    'timestamp','keepthis',0);
if pickingFactor == 20
    XLO = computeQuadratures(data8bitLO(:,:,1:3),configLO, ...
        CALIBRATION_CH1,'LocationOffset',8, ...
        'IntegrationWindow',30,'MinPeakDistance',75,'DutyCycle',dutycycle);
else
    XLO = computeQuadratures(data8bitLO(:,:,1:3),configLO, ...
    CALIBRATION_CH1,'DutyCycle',dutycycle);
end

if strcmp(corrRemove,'yes')
        %Removing correlations with precedent pulses
        XLO(:,:,1:3) = bsxfun(@minus, XLO(:,:,1:3), mean(XLO(:,:,1:3)));
        dispstat('Removing correlations from LO... ','timestamp','keepthis',0);
        XLO = correlationCompensation(XLO);
end

% Calculate the variance piece-wise to compensate slow drifts (e.g.
% piezos)
NLO = mean(var(XLO));

% Compute quadratures for target quantum state
dispstat('Computing target quadratures for Channels 1 to 3', ...
    'timestamp','keepthis',0);
if pickingFactor == 20
    X = computeQuadratures(data8bitSIG(:,:,1:3),configSIG, ...
        CALIBRATION_CH1,'LocationOffset',8, ...
        'IntegrationWindow',30,'MinPeakDistance',75,'DutyCycle',dutycycle);
else
    X = computeQuadratures(data8bitSIG(:,:,1:3),configSIG, ...
        CALIBRATION_CH1,'DutyCycle',dutycycle);
end

for iCh = 1:3
    % Calibration of quadratures to vacuum state
    X(:,:,iCh) = Norm * X(:,:,iCh) / sqrt(NLO(iCh));
    
    % Removing Offsets
    dispstat(['Removing piecewise offset from channel ',num2str(iCh), ...
        ' ...'],'timestamp','keepthis',0);
    X(:,:,iCh) = bsxfun(@minus, X(:,:,iCh), mean(X(:,:,iCh)));
    
    if strcmp(corrRemove,'yes')
        %Removing correlations with precedent pulses
        dispstat(['Removing correlations from channel ',num2str(iCh), ...
            ' ...'],'timestamp','keepthis',0);
        X(:,:,iCh) = correlationCompensation(X(:,:,iCh));
    end
    
    % Cut the raw data into segments of equal length according to piezo
    % modulation
    dispstat(['Reshaping channel ',num2str(iCh), ...
        ' into piezo segments ...'],'timestamp','keepthis',0);
    [Y,piezoSign] = piezoSegments(timestamps,X(:,:,iCh),'cut');
    switch iCh
        case 1
            X1 = Y;
        case 2
            X2 = Y;
        case 3
            X3 = Y;
    end
end % iCh
end % prepare3ChData

