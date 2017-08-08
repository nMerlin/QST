function [ theta, radPerUmUp, radPerUmDown ] = piezoCalibration( data8bit, timestamps, config )
%PIEZO-CALIBRATION computes the phaseshift per micrometer of a piezo motor
%
%   The dataset should be recorded with strong interference of the two
%   beams in a Mach-Zehnder-Interferometer on a photodiode while moving the
%   piezo actuator in a triangular pattern. The data is only recorded while
%   inside a specific range of the triangle given by the TRIGGERBUFFER
%   parameters in the CONFIG file. Then, the calibration values for the
%   upward movement RADPERUMUP and downward movement RADPERUMDOWN will be
%   returned.

SAMPLERATE_HZ = ...
    config.SpectrumCard.Clock.SamplingRate0x28MHz0x29_DBL(1) * 1e6;

PIEZO_DISTANCE_UM = (1 - config.E725.Piezo3.MaxTriggerBuffer_DBL - ...
    config.E725.Piezo3.MinTriggerBuffer_DBL) * ...
    config.E725.Piezo3.Amplitude_DBL;

% Massive averaging of the data8bit-matrix for easier handling and low-pass
% filtering
[nRows,nColumns] = size(data8bit);
modulation = zeros(1,nColumns);
t = zeros(1,nColumns);
deltaT = nRows * 1/SAMPLERATE_HZ * 0.5;
for iColumn=1:nColumns
    modulation(iColumn) = sum(data8bit(:,iColumn))/nRows;
    t(iColumn) = deltaT + double(timestamps(iColumn))/SAMPLERATE_HZ;
end

% Dividing the dataset in segments corresponding to single flanks of the
% triangular modulation
segments = piezoSegments(timestamps, modulation);
tSegments = piezoSegments(timestamps, t);

% Fit and Plot the data segments
nSegments = size(segments,2);
lSegments = size(segments,1);
radPerUm = NaN(1,nSegments);
theta = NaN(size(segments));
hold on;
for iSegment=1:nSegments
    % Fit sine function to the data
    x = tSegments(:,iSegment);
    y = segments(:,iSegment);
    [fitParams, fitFunction] = fitSinusoidal(x, y);
    
    % Phase values (discretized in 64 steps, may rework this)
    thetaColumn = round( ...
        mod(2*pi/fitParams(2)*x + 2*pi/fitParams(3) + pi/2,2*pi),1);
    theta(:,iSegment) = [thetaColumn ...
        NaN(1, lSegments-length(thetaColumn))];
    
    % Phase calibration values
    phasePeriod = fitParams(2);
    piezoVelocity = PIEZO_DISTANCE_UM/(x(length(x(~isnan(x))))-x(1));
    radPerUm(iSegment) = 2*pi/(phasePeriod*piezoVelocity);
    plot(x, y , 'b', x, fitFunction(fitParams,x), 'r');
end
hold off;

% Distinguish upward and downward movement
isFirstSegmentUp = (tSegments(1,2)-tSegments(1,1) > ...
    tSegments(1,3)-tSegments(1,2));
radPerUmUp = mean(radPerUm((2-isFirstSegmentUp):2:nSegments));
radPerUmDown = mean(radPerUm(1+isFirstSegmentUp:2:nSegments));

% Saving important Variables (to delete the raw data manually)
dateString = datestr(datetime('now'),'yyyymmddTHHMMSS');
save(['piezoCalibration-' dateString '.mat'], ...
    'radPerUmUp', 'radPerUmDown', 'radPerUm', 'segments','tSegments', ...
    'config','timestamps');

end

