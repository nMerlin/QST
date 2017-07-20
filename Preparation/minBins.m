function [bins, nBins, edges] = minBins(A,varargin)
%MINBINS Compute bin centers for minimum bin distance in A
%
% Note: The input array A should already exhibit some sort of
% coarse discretization.

%% Optional input paramters
optArgs = {'asym'};
optArgs(1:length(varargin)) = varargin;
[symOpt] = optArgs{:};

uniq = unique(A(:));
hDisc = min(diff(uniq(diff(uniq)>1000*eps)));

%% Compute bin centers and edges
% There are two options: The bins can either be calculated according to the
% maximum and minimum of the data, or symmetrically around 0.
if strcmp(symOpt,'asym')
    bins = min(uniq):hDisc:max(uniq); % bin centers
    edges = (min(uniq)-hDisc/2):hDisc:(max(uniq)+hDisc/2);
elseif strcmp(symOpt,'sym')
    maxVal = max(-min(uniq),max(uniq));
    bins = -maxVal:hDisc:maxVal;
    edges = (-maxVal-hDisc/2):hDisc:(maxVal+hDisc/2);
end
nBins = numel(bins);

end
