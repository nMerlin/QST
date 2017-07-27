function [variance, bin, XOut,meanXBinned,varXBinned ] = assessTheta(theta, X)


%% Create and plot histogram of phase 
% The number of bins is 1000 by default.
xHist = linspace(0,2*pi,1000);
disc = min(diff(xHist));
histEdges = (0-disc/2):disc:(2*pi+disc/2);


h = histogram(theta,histEdges,'Normalization','pdf');
%plotting
h.EdgeColor = 'b';

axis([0 2*pi 0 max(h.Values)+0.05]);
xlabel('\theta');
ylabel('probability density');

%variance of observations pf phase per bin --> Shows if phase is equally 
% distributed
variance = var(h.Values);

%% Histogram of quadrature values per phase 
nIntervals = 1000;
[N,~,bin] = histcounts(theta,nIntervals);
[~,I] = sort(bin);
X = X(I);

XOut = NaN(ceil(3*length(theta)/nIntervals), nIntervals);
for iInterval = 1 : nIntervals
    start = 1+sum(N(1:iInterval-1));
    stop = start+N(iInterval)-1;
   % thetaOut(1:N(iInterval),iInterval) = theta(start:stop);
    XOut(1:N(iInterval),iInterval) = X(start:stop);
end
meanXBinned = mean(XOut, 'omitnan');
varXBinned = var(XOut, 'omitnan');

end
