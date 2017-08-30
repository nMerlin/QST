function [phaseVariance, varXvar] = assessTheta(theta, X, varargin)
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
Norm = 1/sqrt(2);
%Norm ist the factor in the relation between the quadratures and the ladder
%operators: q = Norm*(a^{+} + a), p = Norm*i*(a^{+} - a)
%typical values are 1/sqrt(2) or 1/2. 


%% Validate and parse input arguments
p = inputParser;
% How many bins will be used to examine the distribution of phase values?
defaultPhaseBins = 180;
% How many bins will be used to examine the variance in quadrature values?
% If VarBins is too high, there are bins without or with only a few X
% Values, which results in a high variance of the X variance.
defaultVarBins = defaultPhaseBins;
addParameter(p,'PhaseBins',defaultPhaseBins,@isnumeric);
addParameter(p,'VarBins',defaultVarBins,@isnumeric);
parse(p,varargin{:});
c = struct2cell(p.Results);
[phaseBins,varBins] = c{:};

%% Strip X and theta of their NaN values
X = X(~isnan(X));
theta = theta(~isnan(theta));

%% Create and plot histogram of phase 
% The number of bins is 1000 by default.
histEdges = linspace(0,2*pi,phaseBins+1);

h = histogram(theta,histEdges,'Normalization','pdf');
% Normalization 'pdf' means 'Probability density function estimate'

%variance of observations of phase per bin --> Shows if phase is equally 
% distributed
phaseVariance = var(h.Values);

h.EdgeColor = 'b';
axis([0 2*pi 0 max(h.Values)+0.05]);
xlabel('\theta');
ylabel('Estimation of Probability Density');
title('Distribution of Phase Values');
text('Units','Normalized','Position',[0.6,0.9],'String', ...
    ['Variance = ' num2str(phaseVariance)],'EdgeColor','k');



%% Sorting of Quadrature Values into Phase Bins 
[N,varEdges,bin] = histcounts(theta,varBins);
[~,I] = sort(bin);
X = X(I);

XOut = NaN(max(N), varBins);
for iInterval = 1 : varBins
    start = 1+sum(N(1:iInterval-1));
    stop = start+N(iInterval)-1;
    XOut(1:N(iInterval),iInterval) = X(start:stop);
end
% Compute mean and variance of Quadrature Values for each phase bin
meanXBinned = mean(XOut, 'omitnan');
varXBinned = var(XOut, 'omitnan');
% Measure how constant the variance is
varXvar = var(varXBinned, 'omitnan');

%waitforbuttonpress;
%clf();
figure
hold on;
xAxis = varEdges(1:end-1)+min(diff(varEdges))/2;
xHist = histEdges(1:end-1)+min(diff(histEdges))/2;
plot(xAxis,meanXBinned,'b.',xAxis,varXBinned,'r-o');
plot(xAxis,Norm^2*ones(varBins),'k.','lineWidth',0.5);
plot(xHist,h.Values*10,'g-o');
legend('Mean of phase-binned X','Variance of phase-binned X', ...
    'Variance for coherent state');
xlabel('\theta');
text('Units','Normalized','Position',[0.3,0.1],'String', ...
    ['Variance of Variance = ' num2str(varXvar)],'EdgeColor','k');
title('Phase-Binned Quadrature Values');
set(gca,'YLim',[min(meanXBinned) max(varXBinned)+2]);

end
