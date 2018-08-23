function [H, binsA, binsB] = histogram2D(A,B,varargin)
%HISTOGRAM2D Compute the 2D histogram of the data given in A and B
%
% Optional input arguments:
%   histogram2D(A,B,'plot'): Compute the histogram and plot its heat map.
%       Default value is 'noplot'.
%   histogram2D(A,B,~,binsA,binsB):
%       If binsA and binsB are scalars: Have exactly binsA bins in the
%       A-direction and binsB bins in the B-direction. Default value is
%       1000. 0 will result in as many bins as are necessary to retain
%       existing discretization (only useful for previously discretized
%       data).
%       If binsA and binsB are vectors: Assume that they describe already
%       the binning.

%% Handle optional input arguments
nVarargin = length(varargin);
optArgs = {'noplot' 1000 1000};
optArgs(1:nVarargin) = varargin;
[plotArg,nBinsA,nBinsB] = optArgs{:};

%% Calculate binning and histogram
A = A(:); B = B(:);
if length(nBinsA)==1 && length(nBinsB)==1
    if nBinsA == 0
        [binsA, nBinsA] = minBins(A,'sym');
    else
        maxVal = max(-min(A),max(A));
        binsA = linspace(-maxVal,maxVal,nBinsA);
    end
    if nBinsB == 0
        [binsB, nBinsB] = minBins(B,'sym');
    else
        maxVal = max(-min(A),max(A));
        binsB = linspace(-maxVal,maxVal,nBinsB);
    end
elseif length(nBinsA)>1 && length(nBinsB)>1
    binsA = nBinsA;
    binsB = nBinsB;
    nBinsA = length(binsA);
    nBinsB = length(binsB);
else
    error('Mixed input data types for nBinsA and nBinsB!');
end
% map A and B to bin indices
Ai = round(interp1(binsA,1:nBinsA,A,'linear','extrap'));
Bi = round(interp1(binsB,1:nBinsB,B,'linear','extrap'));
% limit indices to [1,nBins]
Ai = max(min(Ai,nBinsA),1);
Bi = max(min(Bi,nBinsB),1);
% Count number of elements in each bin
H = accumarray([Ai(:) Bi(:)], 1, [nBinsA nBinsB]);

%% Plot
if strcmp(plotArg,'plot')
    imagesc(binsA,binsB,H), axis on
    colormap hot; colorbar;
end

end