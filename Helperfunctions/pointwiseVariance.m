function [ locs, pvar ] = pointwiseVariance( data, varargin )
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
%
%   [___] = POINTWISEVARIANCE(___,'showplot') additionally plots the signal
%   and overlays the peak values

%% Handle optional input arguments
verbose = 0;
showplot = 0;
quiet = 'notquiet';
if nargin > 1
    for i = 1:nargin-1
        eval([varargin{i} '=1;']);
    end
end
if verbose == 0
    quiet = 'quiet';
end

%% Subtract DC offset and calculate pointwise variance
assert(ismatrix(data),'DATA is not a matrix!');
[rows,columns] = size(data);

for column=1:columns
    data(:,column) = data(:,column) - mean(data(:,column));
end

pvar = zeros(rows,1);
for row=1:rows
    pvar(row) = var(double(data(row,:)));
end
% The following vectorized calculation needs 16 GB more memory than the
% preceding for-loop, if the int8 input matrix 'data' has a size of 1 GB.
% This is because the type-casting 'double(data)' creates a new matrix with
% size 8GB (int8->double) and the 'var' command needs to make a copy of
% this matrix. However, the for loop is significantly slower than the
% vectorized calculation ('prepare3ChData' takes 1min 44s in the vectorized
% version vs. 3min 17 with the for-loop.
%pvar = var(double(data),0,2);

%% Find peaks
options.MinPeakHeight = mean(pvar);
options.MinPeakDistance = 10;
[~,locs] = findpeaks(pvar,options);
if showplot == 1
    findpeaks(pvar,options);
end
assert(not(isempty(locs)),'No maxima found in the pointwise variance!');

end

