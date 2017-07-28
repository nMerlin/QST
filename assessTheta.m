function [phaseVariance, varXvar ] = assessTheta(theta, X)
%This function assesses the final result of X3 and theta. It plots the
%distribution of the computed phase and the variance of the distribution.
%It also sorts the X values into phase bins.

%Input parameters:
% theta, X - selected phase and quadrature values obtained with
% selectRegion. 
%
%Output parameters:
%- phaseVariance: variance of the histogram values of the phase, i. e. of
% the probability density.
%- varXvar: The variance of the variance of the X values in each phase bin. 

%% Create and plot histogram of phase 
% The number of bins is 1000 by default.
xHist = linspace(0,2*pi,1000);
disc = min(diff(xHist));
histEdges = (0-disc/2):disc:(2*pi+disc/2);

h = histogram(theta,histEdges,'Normalization','pdf');

%variance of observations of phase per bin --> Shows if phase is equally 
% distributed
phaseVariance = var(h.Values);

h.EdgeColor = 'b';
axis([0 2*pi 0 max(h.Values)+0.05]);
xlabel('\theta');
ylabel('probability density');
title('Phase Distribution');
text('Units','Normalized','Position',[0.6,0.9],'String',['Variance = ' num2str(phaseVariance)],'EdgeColor','k');



%% Sorting of Quadrature Values into Phase bins 
nIntervals = 1000;
[N,~,bin] = histcounts(theta,nIntervals);
[~,I] = sort(bin);
X = X(I);

XOut = NaN(ceil(3*length(theta)/nIntervals), nIntervals);
for iInterval = 1 : nIntervals
    start = 1+sum(N(1:iInterval-1));
    stop = start+N(iInterval)-1;
    XOut(1:N(iInterval),iInterval) = X(start:stop);
end
%compute mean and variance of Quadrature Values for each phase bin
meanXBinned = mean(XOut, 'omitnan');
varXBinned = var(XOut, 'omitnan');
%Measure how constant the variance is
varXvar = var(varXBinned, 'omitnan');

waitforbuttonpress;
clf();
plot(1:nIntervals,meanXBinned,'b.',1:nIntervals,varXBinned,'r.');
legend('Mean of phase-binned X','Variance of phase-binned X');
xlabel('Number of Bin');
text('Units','Normalized','Position',[0.1,0.1],'String',['Variance of Variance = ' num2str(varXvar)],'EdgeColor','k');
title('Phase-Binned Quadrature Values');
end
