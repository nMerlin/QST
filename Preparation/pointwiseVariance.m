function [ locs, pvar ] = pointwiseVariance( data )
% POINTWISEVARIANCE  Calculates the pointwise variance of the given data
% arrays. The input DATA is a matrix and each column is treated as a
% separate data array. The algorithm removes the dc-offset before the
% calculation.
%
%   LOCS = POINTWISEVARIANCE(DATA) Returns the locations of the pointwise
%   variance maxima (as indices) that sit above the mean value.
%
%   [LOCS,PVAR] = POINTWISEVARIANCE(DATA) Additionally outputs the array
%   PVAR which consists of the pointwise variance data points.

assert(ismatrix(data),'DATA is not a matrix!');

[rows,columns] = size(data);
pvar = zeros(rows,1);

for column=1:columns
    data(:,column) = data(:,column) - mean(data(:,column));
end

for row=1:rows
    pvar(row) = var(double(data(row,:)));
end

[~,locs] = findpeaks(pvar,'MinPeakHeight',max(pvar)*0.15,'MinPeakDistance',40);
assert(not(isempty(locs)),'No maxima found in the pointwise variance!');

end

