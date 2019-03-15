function X = preparePhAvData(filenameLO,filenameSIG,varargin)
%PREPAREPHAVDATA Process data from phase-averaged 1-channel measurements
%
% Optional Input Arguments:
%   'CorrRemove': When 'yes', correlations from preceding pulses will be
%     removed.

%% Validate and parse input arguments
p = inputParser;
defaultCorrRemove = 'no';
defaultChannel = 1;
addParameter(p,'CorrRemove',defaultCorrRemove);
addParameter(p,'Channel',defaultChannel,@isnumeric);
parse(p,varargin{:});
c = struct2cell(p.Results);
[channel,corrRemove] = c{:};

CALIBRATION_CH1 = 4.596047840078126e-05; % Ampere per Volt

%Norm ist the factor in the relation between the quadratures and the ladder
%operators: q = Norm*(a^{+} + a), p = Norm*i*(a^{+} - a)
%typical values are 1/sqrt(2) or 1/2.
Norm = 1/sqrt(2);

% Compute number of LO photons
dispstat('Computing number of LO photons ...','timestamp','keepthis',0);
[data8bit,config,~]=load8BitBinary(filenameLO,'dontsave');
XLO = computeQuadratures(data8bit(:,:,channel), config, CALIBRATION_CH1);

if strcmp(corrRemove,'yes')
    %Removing correlations with preceding pulses
    XLO = bsxfun(@minus, XLO, mean(XLO));
    dispstat('Removing correlations from LO ...','timestamp','keepthis',0);
    XLO = correlationCompensation(XLO);
end

% Calculate the variance piece-wise to compensate slow drifts (e.g. piezos)
NLO = mean(var(XLO));

% Compute quadratures for target quantum state
dispstat('Computing quadratures for target quantum state ...',...
    'timestamp','keepthis',0);
[data8bit,config,~]=load8BitBinary(filenameSIG,'dontsave');
X = computeQuadratures(data8bit(:,:,channel), config, CALIBRATION_CH1);

if strcmp(corrRemove,'yes')
    %Removing correlations with preceding pulses
    X = bsxfun(@minus, X, mean(X));
    dispstat('Removing correlations from signal ...', ...
        'timestamp','keepthis',0);
    X = correlationCompensation(X);
end

% Calibration of quadratures to vacuum state
X = Norm * X / sqrt(NLO);

end

