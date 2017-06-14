function [nAv, width, dip, distTherm, distCoh] = plotPhAvData(X, filename, varargin)
%PLOTPHAVDATA Evaluates the quadratures X from a phase-averaged state
%   
%   Mandatory Paramters:
%       X - Matrix containing the measured quadrature values. You can
%           create it for example with the function PREPAREPHAVDATA
%       FILENAME - this string will be included in the filename of 
%           the png-file

% Optional input arguments
verbose = 0;
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

% Remove offset
X = X(:);
X = X - mean(mean(X));

% Mean photon number
nAv = mean(X.^2)-0.5;

% Adjust histogram edges to discretization and plot histogram
uniq = unique(X(:));
maxValue = max(-min(uniq),max(uniq));
hDisc = min(diff(uniq)); % discretization
histEdges = (-maxValue-hDisc/2):hDisc:(maxValue+hDisc/2);
xHist = min(uniq):hDisc:max(uniq);
h = histogram(X,histEdges,'Normalization','probability');
hold on;
h.EdgeColor = 'b';

%%% Compute and plot theory curves
% Thermal State
xAxis = -maxValue:hDisc:maxValue;
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
maxX = 2*sqrt(nAv)+3;
xlim([-maxX maxX]);

% Calculate dip height
dip = abs(max(h.Values)-h.Values(floor(length(h.Values))));

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

