function [uniformX,uniformTheta] = seriesUniformSampling(selX,selTheta,varargin)
%SERIESUNIFORMSAMPLING Uniformly sample in 'selTheta'.

%% Validate and parse input arguments
p = inputParser;
defaultNBins = 100;
addParameter(p,'NBins',defaultNBins,@isnumeric);
parse(p,varargin{:});
c = struct2cell(p.Results);
[nBins] = c{:};

[thetaN,~,bin] = histcounts(selTheta,nBins);
nMin = min(thetaN);
[uniformX,uniformTheta] = deal(zeros(nBins*nMin,1));
for iBin = 1:nBins
    r = randi(thetaN(iBin),nMin,1);
    X = selX(bin==iBin);
    Theta = selTheta(bin==iBin);
    uniformX(((iBin-1)*nMin+1):iBin*nMin) = X(r);
    uniformTheta(((iBin-1)*nMin+1):iBin*nMin) = Theta(r);
end

end

