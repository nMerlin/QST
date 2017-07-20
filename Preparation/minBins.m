function [bins, nBins, edges] = minBins(A)
%MINBINS Compute bin centers for minimum bin distance in A
%
% Note: The input array A should already exhibit some sort of
% coarse discretization.

uniq = unique(A(:));
hDisc = min(diff(uniq(diff(uniq)>1000*eps)));
bins = min(uniq):hDisc:max(uniq); % bin centers
nBins = numel(bins);
edges = (min(uniq)-hDisc/2):hDisc:(max(uniq)+hDisc/2);

end
