function [ X, theta ] = mainPrepareData( filename )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

CALIBRATION_CH1 = 3.462042300029793e+03;
NLO = 4.379815002805815e+08;

[data8bit,config,timestamps]=load8BitBinary(filename,'dontsave');
[locs,~] = pointwiseVariance(data8bit);
[~, X] = correlation(0, data8bit, locs);
X = piezoSegments(timestamps, X);

% Calibration of quadratures
X = X * CALIBRATION_CH1;
X = X / sqrt(2 * NLO);

% Compute relative phases
[X, theta] = computeTheta(X);

% Remove unphysical offset
for iSeg = 1 : size(X,2)
    Seg = X(:,iSeg);
    X(:,iSeg) = Seg - (max(Seg) - 0.5 * (max(Seg) - min(Seg)));
end

% Saving important Variables (to delete the raw data manually)
dateString = datestr(datetime('now'),'yyyymmddTHHMMSS');
save(['quadratureDataset-' dateString '.mat'], ...
    'X', 'theta');

end