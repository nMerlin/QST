function [ X, theta ] = mainPrepareData( filenameLO, filenameSIG )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

CALIBRATION_CH1 = 4.596047840078126e-05; % Ampere per Volt

% Compute number of LO photons
[data8bit,config,~]=load8BitBinary(filenameLO,'dontsave');
XLO = computeQuadratures(data8bit, config, CALIBRATION_CH1);
NLO = mean(var(XLO));

% Compute quadratures for target quantum state
[data8bit,config,timestamps]=load8BitBinary(filenameSIG,'dontsave');
X = computeQuadratures(data8bit, config, CALIBRATION_CH1);
X = piezoSegments(timestamps, X);
% Calibration of quadratures
X = X / sqrt(2 * NLO);

% Compute relative phases and removes offset
[X, theta] = computeTheta(X);

% Saving important Variables (to delete the raw data manually)
dateString = datestr(datetime('now'),'yyyymmddTHHMMSS');
save(['quadratureDataset-' dateString '.mat'], ...
    'X', 'theta');

end