function [X, piezoSign] = preparePicked1ChData(filenameLO, filenameSIG, varargin)
%PREPARE1CHDATA Returns quadratures of a 1-Channel-Measurement
%
% The quadrature matrices are already cut into piezo segments.
%
%   filenameLO - filename of the LO-data used for correct normalization
%   filenameSIG - filename of the raw 3-Channel-Data
%
% Output Arguments:
%   piezoSign: +1 means piezo moves in the first segment from 0 to 2 um
%              -1 means piezo moves in the first segment from 2 to 0 um

%% Validate and parse input arguments
p = inputParser;
defaultChannel = 1;
addParameter(p,'Channel',defaultChannel,@isnumeric);
parse(p,varargin{:});
c = struct2cell(p.Results);
[channel] = c{:};

CALIBRATION_CH1 = 4.596047840078126e-05; % Ampere per Volt

%Norm ist the factor in the relation between the quadratures and the ladder
%operators: q = Norm*(a^{+} + a), p = Norm*i*(a^{+} - a)
%typical values are 1/sqrt(2) or 1/2.
Norm = 1/sqrt(2);

dispstat('','init','timestamp','keepthis',0);
dispstat('Loading LO data ...','timestamp','keepthis',0);
[data8bitLO,configLO,~]=load8BitBinary(filenameLO,'dontsave');

dispstat('Loading Signal data ...','timestamp','keepthis',0);
[data8bitSIG,configSIG,timestamps]=load8BitBinary(filenameSIG,'dontsave');

% Compute number of LO photons
dispstat('Computing number of LO photons ...', ...
    'timestamp','keepthis',0);
XLO = computeQuadratures(data8bitLO(:,:,channel), configLO, ...
    CALIBRATION_CH1);

% Calculate the variance piece-wise to compensate slow drifts (e.g.
% piezos)
NLO = mean(var(XLO));

% Compute quadratures for target quantum state
dispstat('Computing target quadratures ...', ...
    'timestamp','keepthis',0);
X = computeQuadratures(data8bitSIG(:,:,channel), configSIG, ...
    CALIBRATION_CH1);

% Calibration of quadratures to vacuum state
X = Norm * X / sqrt(NLO);

% Removing Offset
dispstat('Removing global offset ...','timestamp','keepthis',0);
X = X - mean(mean(X));

% Cut the raw data into segments of equal length according to piezo
% modulation
dispstat('Reshaping quadratures into piezo segments ...', ...
    'timestamp','keepthis',0);
[X,piezoSign] = piezoSegments(timestamps,X,'cut');

end % prepare1ChData

