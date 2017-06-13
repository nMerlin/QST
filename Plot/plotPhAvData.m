function [nAv, width, dip, distTherm, distCoh] = plotPhAvData(X, filename)
%PLOTPHAVDATA Evaluates the quadratures X from a phase-averaged state
%   
%   Input Paramters:
%   X - Matrix containing the measured quadrature values. You can create it
%       for example with the function PREPAREPHAVDATA
%   FILENAME - this string will be included in the filename of the png-file

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
WF = thermWigner(xAxis,xAxis,nAv);
qThermal = sum(WF);
WF = cohWigner(xAxis,xAxis,nAv);
qCoherent = phaseAveragedDist(WF,xAxis);
plot(xAxis,qThermal,'linewidth',2);
plot(xAxis,qCoherent,'linewidth',2);
legend('Measured','Thermal','Coherent');
    
% Calculate euclidean distance from both theoretical states
% start = (length(h.Values) - length(qThermal))/2+1;
% stop = start + length(qThermal) - 1;
distTherm = sqrt(sum((qThermal - h.Values).^2));
distCoh = sqrt(sum((qCoherent - h.Values).^2));

xlabel('Q');
ylabel('Probability');
title({filename,['Phase-averaged quantum state (n=',num2str(nAv),')']});

% Adjust axis limits
width = fwhm(xHist,h.Values);
maxX = 2*sqrt(nAv)+3;
xlim([-maxX maxX]);
hold off;

% Calculate dip height
dip = abs(max(h.Values)-h.Values(floor(length(h.Values))));

% Save figure
print(strcat('PHAV-',filename,'-',num2str(nAv),'photons.png'), '-dpng');

end

