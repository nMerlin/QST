function [ X, theta ] = mainPrepareData( filenameLO, filenameSIG)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Norm = 1/sqrt(2);
%Norm ist the factor in the relation between the quadratures and the ladder
%operators: q = Norm*(a^{+} + a), p = Norm*i*(a^{+} - a)
%typical values are 1/sqrt(2) or 1/2. 

CALIBRATION_CH1 = 4.596047840078126e-05; % Ampere per Volt

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
[X, piezoSign] = piezoSegments(timestamps, X);
% Calibration of quadratures
X = Norm * X / sqrt(NLO);

% Compute relative phases and removes offset
[X, theta] = computeTheta(X,piezoSign,'verbose');

% Saving important Variables (to delete the raw data manually)
dispstat('Saving ...','timestamp','keepthis',0);
dateString = datestr(datetime('now'),'yyyymmddTHHMMSS');
save(['quadratureDataset-' strrep(num2str(filenameSIG),'LOwithLO.raw','') dateString '.mat'], ...
    'X', 'theta');

end