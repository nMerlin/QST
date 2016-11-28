function [ X ] = computeQuadratures( data8bit, config, amperePerVolt )
%COMPUTEQUADRATURES Compute quadrature values in number of photons
%
% CALIBRATION is given in A/V

INTEGRATION_DUTY_CYCLE = 1/3;
SAMPLERATE = config.SpectrumCard.Clock.SamplingRate0x28MHz0x29_DBL * 10e6;
ELEMENTARY_CHARGE = 1.6021766208e-19;

switch config.SpectrumCard.Channel00.Range_I32
    case 0
        INT8_TO_VOLTAGE = 0.200/128;
    case 1
        INT8_TO_VOLTAGE = 0.500/128;
    case 2
        INT8_TO_VOLTAGE = 1.0/128;
    case 3
        INT8_TO_VOLTAGE = 2.5/128;
end

% Identify integration centers
[locs,~] = pointwiseVariance(data8bit);

% Eliminate locations whose corresponding window would be outside the range
% of DATA.
[nRows, nColumns] = size(data8bit);
window = round(INTEGRATION_DUTY_CYCLE * mean(diff(locs)));
if (locs(1)<=ceil(window/2))
    locs = locs(2:end);
end
if ((nRows-locs(end))<ceil(window/2))
    locs = locs(1:length(locs)-1);
end

% Integration loop
nWindows = length(locs);
X = zeros(nWindows, nColumns);
for iColumn = 1 : nColumns
    for iWindow = 1 : nWindows
        % Integration and calibration step
        start = locs(iWindow)-ceil(window/2);
        stop = locs(iWindow)+ceil(window/2);
        windowTime = (stop-start+1) * 1 / SAMPLERATE;
        X(iWindow, iColumn) = sum(data8bit(start:stop, iColumn)) * ...
            INT8_TO_VOLTAGE * amperePerVolt * windowTime / ...
            ELEMENTARY_CHARGE;
    end % iWindow
end % iColumn

end

