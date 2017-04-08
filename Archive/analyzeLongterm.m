function [t, thetaUp, thetaDown] = analyzeLongterm(data8bit, timestamps)
%ANALYZELONGTERM Compute Phase/Time plot for long measurements (e.g. 300s)
%
%   [t, thetaUp, thetaDown] = analyzeLongterm(data8bit, timestamps)
%   
%   timestamps: timestamps of double-trigger events (piezo trigger HIGH and
%               main LO trigger POS FLANK)
%   data8bit: matrix (nSamples x nRecords x nChannels)
%   t: time in samples
%   thetaUp: Computed phase values on upward piezo flank
%   thetaDown: Computed phase values on downward piezo flank
%   CAUTION: Directions of Up and Down could be reversed
%
%   Function written for dataset '2017-03-28-Phase-Noise-Longterm' under
%   the following measurement conditions:
%       - Only CH1 and CH2 on ADC used with sampling rate 305175 Hz
%       - CH1 was connected to HD-Det. RevC SN02 and CH2 to HD-Det. RevC
%         SN03.
%       - Both LO powers were 5mW
%       - Both SIG originated from the same laser as LO
%       - Both piezo actuators were driven at 50Hz with 2um amplitude in a
%         triangular fashion and the piezo controller sent a trigger signal
%         to the secondary trigger input of the ADC.
%       - The ADCs main trigger was connected to a Photodiode monitoring
%         the LO.

[nSamples, nRecords, ~] = size(data8bit);

t = zeros(nSamples, nRecords, 'uint64');

for iSegments = 1:nRecords
    t(:,iSegments) = uint64(0:95) + timestamps(iSegments);
end

% Analyze only selected channel
channel = 1;
data = data8bit(:,:,channel);

% Dividing into PiezoSegments
data = piezoSegments(timestamps, data);
t = piezoSegments(timestamps, t);

[~,~,nSegments] = size(data);

% Averaging over 16 consecutive samples
nAverage = 16;
data = mean(reshape(data,nAverage,[],nSegments));
t = mean(reshape(t,nAverage,[],nSegments));

data = reshape(data,[],nSegments);
t = reshape(t,[],nSegments);
t = mean(t(:,1:2:end));

[~, nSegments] = size(data);
theta = zeros(nSegments,1);
fvals = zeros(nSegments,1);
exitflags = zeros(nSegments,1);

% Reduce computational load by reducing number of fits
%nSegments = 100;

tic

parfor iSeg = 1:nSegments
%for iSeg = nSegments:nSegments
    yFit = data(16:end-32,iSeg);
    xFit = 0:length(yFit)-1;
    
    xFit(isnan(yFit)) = NaN;

    % Somehow the fit only works for the correct x magnitude
    xFitMagnitude = ceil(log10(max(xFit)));
    xFit = xFit / 10^(xFitMagnitude); % scale x-axis for fitting routine
    
    % Fitting
    [fitParams,fvals(iSeg),exitflags(iSeg)] = fitSinusoidal(xFit, yFit, 'rmLin');
    
    % Calculate phase values from fit results
    theta(iSeg) = mod(2 * pi / fitParams(3) + pi / 2, 2*pi);
end
toc

% Clean theta: Remove phase values where the fit didn't converge or the
% fitting error is 10x bigger than average.
meanfval = mean(fvals);
for iSeg = 1:nSegments
    if (exitflags(iSeg) ~= 1) || (fvals(iSeg) > 10*meanfval)
        theta(iSeg) = 0;
    end
end

% Update theta for continuity: Split in the two piezo movement directions
thetaUp = smoothTheta(theta(1:2:end));
thetaDown = smoothTheta(theta(2:2:end));

plotLongterm(t, thetaUp, thetaDown);

end

function [ theta ] = smoothTheta(theta)
    %SMOOTHTHETA converts phase values from [0,2pi] into a continuous
    %trace
    for iTheta = 2:length(theta)
        lastTheta = theta(iTheta-1);
        if theta(iTheta) == 0
            theta(iTheta) = lastTheta;
        else
            nShift = round((lastTheta - theta(iTheta))/(2*pi));
            theta(iTheta) = theta(iTheta) + nShift * 2 * pi;
        end
    end
end

function plotLongterm( time, thetaUp, thetaDown )
%PLOTLONGTERM plot phase values against time

figure('units','normalized','position',[.1 .1 .6 .4]);
plot(time, thetaUp/pi, time, thetaDown/pi);
ylabel('Phase [\pi]');
xlabel('Time [s]');
title('02-Ch23.raw | Ch2 - Ch3');
legend('Uphill','Downhill','Location','Northwest');

end