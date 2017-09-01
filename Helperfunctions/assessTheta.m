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
defaultPhaseBins = 100;
% How many bins will be used to examine the variance in quadrature values?
% If VarBins is too high, there are bins without or with only a few X
% Values, which results in a high variance of the X variance.
defaultVarBins = 100;
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

fig = figure('Units','centimeters','Position',[1,1,14.8,21]);
ax1 = subplot(2,1,1);
h = histogram(ax1,theta,histEdges); % 'Normalization','pdf');
% Normalization 'pdf' means 'Probability density function estimate'

%variance of observations of phase per bin --> Shows if phase is equally 
% distributed
phaseVariance = var(h.Values);
Nbar = length(theta)/phaseBins;

h.EdgeColor = 'b';
axis([0 2*pi 0 max(h.Values)+0.05]);
xlabel('\theta');
ylabel('N(\theta)');
title(['Distribution of Phase Values (',num2str(phaseBins),' Bins)']);
legend(['mean(N(\theta)) = ',num2str(Nbar),'; Var(N(\theta))) = ', ...
    num2str(phaseVariance)]);
%text('Units','Normalized','Position',[0.6,0.9],'String', ...
%    ['Variance = ' num2str(phaseVariance)],'EdgeColor','k');

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
meanVarX = mean(varXBinned);
% Measure how constant the variance is
varXvar = var(varXBinned, 'omitnan');

%waitforbuttonpress;
%clf();
%figure
%hold on;
ax2 = subplot(2,1,2);
hold on;
xAxis = varEdges(1:end-1)+min(diff(varEdges))/2;
plot(xAxis,meanXBinned,'b.',xAxis,varXBinned,'r-o');
plot(xAxis,Norm^2*ones(varBins),'k.','lineWidth',0.5);
legend('Mean of phase-binned X values', ...
    ['mean(Var(X)) = ',num2str(meanVarX),', Var(Var(X)) = ',num2str(varXvar)], ...
    ['Var(Coherent State ) = ',num2str(Norm^2)]);
xlabel('\theta');
title(['Phase-Binned Quadrature Values (',num2str(varBins),' Bins)']);
set(gca,'YLim',[min(meanXBinned) max(varXBinned)+3],'XLim', ...
    [min(xAxis) max(xAxis)]);
fig.PaperPositionMode = 'auto';
print(['Variance-',num2str(varBins),'-varBins-',num2str(phaseBins), ...
    '-phaseBins','.pdf'],'-dpdf');

end
