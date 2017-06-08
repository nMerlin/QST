function [ X ] = preparePhAvData( filenameLO, filenameSIG )
%PREPAREPHAVDATA Process data from phase-averaged 1-channel measurements
%   Detailed explanation goes here

CALIBRATION_CH1 = 4.596047840078126e-05; % Ampere per Volt

%Norm ist the factor in the relation between the quadratures and the ladder
%operators: q = Norm*(a^{+} + a), p = Norm*i*(a^{+} - a)
%typical values are 1/sqrt(2) or 1/2.
Norm = 1/sqrt(2);

% Compute number of LO photons
dispstat('Computing number of LO photons ...','timestamp','keepthis',0);
[data8bit,config,~]=load8BitBinary(filenameLO,'dontsave');
XLO = computeQuadratures(data8bit, config, CALIBRATION_CH1);

% Calculate the variance piece-wise to compensate slow drifts (e.g. piezos)
NLO = mean(var(XLO));

% Compute quadratures for target quantum state
dispstat('Computing quadratures for target quantum state ...',...
    'timestamp','keepthis',0);
[data8bit,config,timestamps]=load8BitBinary(filenameSIG,'dontsave');
X = computeQuadratures(data8bit, config, CALIBRATION_CH1);

% Calibration of quadratures to vacuum state
X = Norm * X / sqrt(NLO);

end

