function [bins, nBins] = minBins(A)
%MINBINS Compute bin centers for minimum bin distance in A
%
% Note: The input array A should already exhibit some sort of
% discretization

uniq = unique(A(:));
hDisc = min(diff(uniq));
bins = min(uniq):hDisc:max(uniq); % bin centers
nBins = numel(bins);

end

