function [nAv, width, dip, distTherm, distCoh] = plotPhAvData(X, filename, varargin)
%PLOTPHAVDATA Evaluates the quadratures X from a phase-averaged state
%   
% Mandatory Paramters:
%   X - Matrix containing the measured quadrature values. You can
%       create it for example with the function PREPAREPHAVDATA
%   FILENAME - this string will be included in the filename of 
%       the png-file
%
% Optional Parameters:
%   'limBins': Limit the number of bins in the histogram to 1000

%% Optional input arguments
verbose = 0;
offset = 0;
limBins = 0;
quiet = 'notquiet';
incoherent = 0;
if nargin > 2
    for i = 3:nargin
        eval([varargin{i-2} '=1;']);
    end
end
if verbose == 0
    quiet = 'quiet';
end

%% Account for different sizes of X
% X can have the shape [nPoints,nRecords] or [nPoints,nRecords,nSegments],
% where _nSegments_ is the number of piezo segments. The following code
% identifies the shape and transforms it to the first type, if necessary.
if ndims(X) == 3
    [nPoints,nRecords,nSegments] = size(X);
    X = reshape(X,[nPoints,nRecords*nSegments]);
end

%%% Piecewise photon number <a^+ a>
if offset == 0
    X = X - mean(mean(X));
end
ada = mean(X.^2)-0.5;

%%% Piecewise <a^+ a^+ a a>, in the following adadaa, and g2(t,t)
adadaa = 2/3*mean(X.^4)-2*ada-0.5;
g2 = adadaa./ada.^2;

%%% Global photon number
% Remove offset
X = X(:);
% Mean photon number
nAv = mean(X.^2)-0.5;

%% Create and plot histogram
% Adjust histogram edges to the discretization of X, if the option
% 'minDisc' is set. Otherwise, the number of bins is 1000 by default.
if limBins == 0
    [xHist,~,histEdges] = minBins(X,'sym');
else
    maxVal = max(-min(X),max(X));
    xHist = linspace(-maxVal,maxVal,1000);
    disc = min(diff(xHist));
    histEdges = (-maxVal-disc/2):disc:(maxVal+disc/2);
end
clf;
h = histogram(X,histEdges,'Normalization','probability');
hold on;
h.EdgeColor = 'b';

%%% Compute and plot theory curves
% Thermal State
xAxis = xHist;
WF_therm = thermWigner(xAxis,xAxis,nAv);
qThermal = sum(WF_therm);
WF_coh = cohWigner(xAxis,xAxis,nAv);
qCoherent = phaseAveragedDist(WF_coh,xAxis);
plot(xAxis,qThermal,'linewidth',2);
plot(xAxis,qCoherent,'linewidth',2);
legend('Measured','Thermal','Coherent');
    
% Calculate euclidean distance from both theoretical states
distTherm = sqrt(sum((qThermal - h.Values).^2));
distCoh = sqrt(sum((qCoherent - h.Values).^2));

xlabel('Q');
ylabel('Probability');
title({filename,['Phase-averaged quantum state (n=',num2str(nAv),')']});

% Adjust axis limits
width = fwhm(xHist,h.Values);
maxX = 3.5*sqrt(nAv)+3;
xlim([-maxX maxX]);

% Calculate dip height
dip = abs(max(h.Values)-h.Values(floor(length(h.Values))));

% Plot inset with piecewise photon numbers
mainFigure = gcf;
insetAx = axes('Parent',mainFigure,'Position',[0.18 0.6 0.15 0.25]);
plot(ada);
set(insetAx,'FontSize',8,'XLim',[1 length(ada)], ...
    'YLim',[min(ada) max(ada)],'XTickLabel','');
title('Photon Number');
insetAx2 = axes('Parent',mainFigure,'Position',[0.18 0.30 0.15 0.25]);
plot(insetAx2, g2);
set(insetAx2,'FontSize',8,'XLim',[1 length(g2)], ...
    'YLim',[min(g2) max(g2)],'XTickLabel','');
xlabel('Time');
title('g2-correlation');
figure(mainFigure);

% Fit incoherent superposition
if incoherent == 1
    minAlpha = 0;
    minDist = distTherm;
    for alpha=0:0.01:1
        WF = alpha*WF_coh+(1-alpha)*WF_therm;
        qMix = phaseAveragedDist(WF,xAxis);
        distMix = sqrt(sum((qMix - h.Values).^2));
        if distMix < minDist
            minDist = distMix;
            minAlpha = alpha;
        end
    end
    WF = minAlpha*WF_coh+(1-minAlpha)*WF_therm;
    qMix = phaseAveragedDist(WF,xAxis);
    plot(xAxis,qMix,'linewidth',2);
end
hold off;

% Save figure
print(strcat('PHAV-',filename,'-',num2str(nAv),'photons.png'), '-dpng');

end

