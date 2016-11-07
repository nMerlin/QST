function [ X, theta ] = mainPrepareData( filename )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

WINDOWSIZE = 40;

[data8bit,config,timestamps]=load8BitBinary(filename,'dontsave');
theta = piezoCalibration(data8bit, timestamps, config);
[locs,~] = pointwiseVariance(data8bit);
[~, X] = correlation(0, data8bit, locs, WINDOWSIZE);
X = piezoSegments(timestamps, X);

% Saving important Variables (to delete the raw data manually)
dateString = datestr(datetime('now'),'yyyymmddTHHMMSS');
save(['quadratureDataset-' dateString '.mat'], ...
    'X', 'theta');

end

