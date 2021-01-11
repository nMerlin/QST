function [uniformX,uniformTheta,uniformIndices] = seriesUniformSamplingIndex(selX,selTheta,varargin)
%SERIESUNIFORMSAMPLING Uniformly sample in 'selTheta'.
%it provides also the indices of the sampled values from the original
%arrays. 

%% Validate and parse input arguments
p = inputParser;
defaultNBins = 100;
addParameter(p,'NBins',defaultNBins,@isnumeric);
parse(p,varargin{:});
c = struct2cell(p.Results);
[nBins] = c{:};

[thetaN,~,bin] = histcounts(selTheta,nBins);
nMin = min(thetaN);
[uniformX,uniformTheta,uniformIndices] = deal(zeros(nBins*nMin,1));
[~,binI] = sort(bin);
stops = cumsum(thetaN);
start = 1;
for iBin = 1:nBins
    r = randi([start,stops(iBin)],nMin,1);
    uniformX(((iBin-1)*nMin+1):iBin*nMin) = selX(binI(r));
    uniformIndices(((iBin-1)*nMin+1):iBin*nMin) = binI(r);
    if nargout>1
        uniformTheta(((iBin-1)*nMin+1):iBin*nMin) = selTheta(binI(r));
    end
    start = stops(iBin)+1;
end

end
