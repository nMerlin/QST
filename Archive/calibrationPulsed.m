function [ calibration ] = calibrationPulsed( vIn, tIn, vOut, tOut, intRange )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

WINDOW = 50e-9;
SAMPLINGRATE_IN = 40e9;
SAMPLINGRATE_OUT = 5e9;
INPUT_RESISTANCE = 100020; % Ohm

% Get WINDOW around one peak in vIn
[~, locs] = findpeaks(vIn,'Annotate', 'extents', ...
    'MinPeakProminence', 0.05, 'MinPeakDistance', 1e-7);
assert(length(locs)==1, 'More than one peak found in vIn!');
targetRange = (1 : (WINDOW * SAMPLINGRATE_IN));
targetRange = targetRange + locs - round(length(targetRange) / 2);
vIn = vIn(targetRange);
tIn = tIn(targetRange) - tIn(locs);

% Get WINDOW around one peak in vOut
if abs(min(vOut)) > abs(max(vOut))
    vOut = -vOut;
end
[~, locs] = findpeaks(vOut,'Annotate', 'extents', ...
    'MinPeakProminence', 0.05, 'MinPeakDistance', 50);
assert(length(locs)==3, 'Not exactly three peaks found in vOut!');
targetRange = (1 : (WINDOW * SAMPLINGRATE_OUT));
targetRange = targetRange + locs(1) - round(length(targetRange) / 2);
vOut = vOut - mean(vOut(tOut < (tOut(locs(1)) - WINDOW / 2)));
vOut = vOut(targetRange);
tOut = tOut(targetRange) - tOut(locs(1));

% Resample vOut
vOut = resample(vOut, SAMPLINGRATE_IN, SAMPLINGRATE_OUT);
tOut = (1 : length(vOut)) * (tOut(end) - tOut(1)) / length(vOut) + tOut(1);

% Optimize overlap
productMax = 0;
iMax = -500;
for iShift = -500 : 500
    product = sum(abs(vOut .* circshift(vIn, iShift)));
    if product > productMax
        iMax = iShift;
        productMax = product;
    end
end % iShift
vIn = circshift(vIn, iMax);

% Compute calibration value
sumOfInputCurrent = 1 / INPUT_RESISTANCE * sum(vIn(intRange));
sumOfOutputVoltage = sum(vOut(intRange));
calibration = sumOfInputCurrent / sumOfOutputVoltage;

end

