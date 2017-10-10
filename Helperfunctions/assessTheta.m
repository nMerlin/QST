function [phaseVariance, varXvar, meanVarX, varTable] = assessTheta(theta, X, varargin)
%This function assesses the final result of X3 and theta. It plots the
%distribution of the computed phase and the variance of the distribution.
%It also sorts the X values into phase bins.

% Input parameters:
%   theta, X - selected phase and quadrature values obtained with
%       selectRegion.
% 
% Optional Input Parameters ('Name','Value'):
%   'Zoom' - Use 'var' to scale the plot of mean quadratures and variance
%       according to the computed variance values. Default: 'none'.
%
%Output parameters:
%- phaseVariance: variance of the histogram values of the phase, i. e. of
% the probability density.
%- varXvar: The variance of the variance of the X values in each phase bin. 
%- meanVarX: The mean of the variance of the X values 
Norm = 1/sqrt(2);
%Norm ist the factor in the relation between the quadratures and the ladder
%operators: q = Norm*(a^{+} + a), p = Norm*i*(a^{+} - a)
%typical values are 1/sqrt(2) or 1/2. 


%% Validate and parse input arguments
p = inputParser;
% How many bins will be used to examine the distribution of phase values?
defaultPhaseBins = 200;
% How many bins will be used to examine the variance in quadrature values?
% If VarBins is too high, there are bins without or with only a few X
% Values, which results in a high variance of the X variance.
defaultVarBins = 200;
defaultHusimi = {};
defaultOutput = 'figure';
defaultZoom = 'none';
defaultFilename = '';
addParameter(p,'PhaseBins',defaultPhaseBins,@isnumeric);
addParameter(p,'VarBins',defaultVarBins,@isnumeric);
addParameter(p,'Husimi',defaultHusimi);
addParameter(p,'Output',defaultOutput);
addParameter(p,'Zoom',defaultZoom,@isstr);
addParameter(p,'Filename',defaultFilename,@isstr);
parse(p,varargin{:});
c = struct2cell(p.Results);
[filename,husimi,output,phaseBins,varBins,zoom] = c{:};

%% Strip X and theta of their NaN values
X = X(~isnan(theta));
theta = theta(~isnan(theta));

%% Create and plot histogram of phase 
% The number of bins is 1000 by default.
histEdges = linspace(0,2*pi,phaseBins+1);

fig = figure('Units','centimeters','Position',[1,1,21,29.7]);
if isempty(husimi)
    ax2 = subplot(2,1,1);
else
    subplot(3,1,1);
    plotHusimi(husimi{1},husimi{2},husimi{3}); % husimi=[O1,O2,iSelect]
    ax2 = subplot(3,1,2);
end
h = histogram(ax2,theta,histEdges); % 'Normalization','pdf');
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
if isempty(husimi)
    ax3 = subplot(2,1,2);
else
    ax3 = subplot(3,1,3);
end
hold on;
xAxis = varEdges(1:end-1)+min(diff(varEdges))/2;
plot(ax3,xAxis,meanXBinned,'b.',xAxis,varXBinned,'r-o');
plot(ax3,xAxis,Norm^2*ones(varBins),'k.','lineWidth',0.5);
legend('Mean of phase-binned X values', ...
    ['Var(X); <Var(X)>=',num2str(meanVarX),', Var(Var(X))=',num2str(varXvar)], ...
    ['Var(Coherent State ) = ',num2str(Norm^2)]);
xlabel('\theta');
title(['Phase-Binned Quadrature Values (',num2str(varBins),' Bins)']);
set(gca,'XLim',[min(xAxis) max(xAxis)]);
if strcmp(zoom,'none')
    set(gca,'YLim',[min(min(meanXBinned),min(varXBinned)), ...
        max(max(varXBinned)+3,max(meanXBinned))]);
elseif strcmp(zoom,'var')
    set(gca,'YLim',[min(varXBinned) max(varXBinned)]);
    %set(gca,'YLim',[0.46 0.54]);
end

%% Output to File
if strcmp(output,'print')
    fig.PaperPositionMode = 'auto';
    print([filename '-Variance-',num2str(varBins),'-varBins-',num2str(phaseBins), ...
        '-phaseBins','.pdf'],'-dpdf');
end

%% Create output table
Phase = xAxis';
meanX = meanXBinned';
varX = varXBinned';
varTable = table(Phase,meanX,varX);

end
