function [bins, nBins] = minBins(A)
%MINBINS Compute bin centers for minimum bin distance in A
%   Detailed explanation goes here

uniq = unique(A(:));
hDisc = min(diff(uniq));
bins = min(uniq):hDisc:max(uniq); % bin centers
nBins = numel(bins);

end

