function [g2vec, ada, meang2, g2Std, XOut, Indices, n] = g2Bins(X, nResolution, varBins, filename,varargin)
% This function first computes the photon numbers from the quadratures X with the
%resolution nResolution. Then it sorts the photon numbers into bins. Then
%it computes g2 for each bin, i.e. for each mean photon number seperately.
%
% Input Arguments:
%   X - Quadratures from a continuous quantum state measurement (equal time
%       spacing between all points)
%   NRESOLUTION - Number of Quadratures used to compute the mean photon
%   number
%   varBins - Amount of Bins into which the photon numbers are sorted
%   range - meang2 will be computed in a range of photon numbers around the
%   middle bin 
%
% Output Arguments:
%   g2vec: g2(0) values sorted according to increasing photon number
%   ada: number of photons, sorted according to increasing photon number
%   meang2: the mean value of g2 in the chosen range of photon numbers
%   g2Std: the standard deviation of the g2 values in the chosen range of photon numbers
%   XOut, Indices, n: Quadratures binned according to photon number and their
%   indices, vector of n for each quadrature

%% Validate and parse input arguments
p = inputParser;
defaultPlotoption = 'yes'; 
addParameter(p,'Plot',defaultPlotoption);
defaultRange = 1/3; 
addParameter(p,'Range',defaultRange, @isnumeric);
parse(p,varargin{:});
c = struct2cell(p.Results);
[plotoption,range] = c{:};

%% Reshaping X according to the starting resolution NRESOLUTION
X = X(:);
nSegments = floor(length(X)/nResolution);
Xres = X(1:nSegments*nResolution);
Xres = reshape(Xres,[nResolution nSegments]); % Consecutive values in columns
X = Xres(:);

%% Piecewise photon number <a^+ a>
Xres = Xres - mean(mean(Xres));
ada = mean(Xres.^2)-0.5;
ada(ada<0) = 0; %set photon numbers < 0 to 0.

%% concatenate them so that we have each element nRes times ... 
n = zeros(size(Xres));
for iRow = 1:nResolution
    n(iRow,1:end) = ada;
end
n = n(:);

%% Make bins of X according to photon number n
[N,~,bin] = histcounts(n,varBins);
[~,I] = sort(bin);
X = X(I);

[XOut,Indices] = deal(NaN(max(N), varBins));
for iInterval = 1 : varBins
    start = 1+sum(N(1:iInterval-1));
    stop = start+N(iInterval)-1;
    XOut(1:N(iInterval),iInterval) = X(start:stop);
    Indices(1:N(iInterval),iInterval) = I(start:stop);
end

%% Compute g2 seperately for each bin:
%Piecewise <a^+ a^+ a a>, in the following adadaa, and g2(t,t)�
ada = mean(XOut.^2,'omitnan')-0.5;
adadaa = 2/3*mean(XOut.^4,'omitnan')-2*ada-0.5;
g2vec = adadaa./ada.^2;
g2vec = g2vec';
g2vec = g2vec(~isnan(g2vec));
ada = ada(~isnan(g2vec));
N=N(~isnan(g2vec));

%% Compute the mean g2 in the middle of the range
if length(g2vec) > 1
    meang2 = mean(g2vec( round(length(g2vec)*(1-range)/2) : round(length(g2vec)*(1+range)/2)));
    g2Std = sqrt(var(g2vec( round(length(g2vec)*(1-range)/2) : round(length(g2vec)*(1+range)/2))));
else
    meang2 = mean(g2vec);
    g2Std = sqrt(var(g2vec));
end

%% plot evaluation
if strcmp(plotoption,'yes')
    subplot(2,1,1);
    plot(ada,g2vec,'o');
    ylim([-0.5 3]);
    xlabel('mean nPhotons');
    ylabel('g2');
    text(0.5, 0.75, ['mean g2 around the middle: ' num2str(meang2,'%.2f')],'Units','normalized');  
    hold on;
    subplot(2,1,2);
    plot(ada, N,'o');
    xlabel('mean nPhotons');
    ylabel('number of values');
    savefig([filename '-g2-Binned-nRes-' num2str(nResolution) '-Bins-'...
        num2str(varBins) '-Range-' num2str(range) '.fig']);
    clf();
end

end

