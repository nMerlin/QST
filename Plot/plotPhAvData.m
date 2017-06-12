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

% Simulations only if nMax<101
nMax = max(ceil(4*nAv),30);
if nMax<101
    % Simulate theoretical thermal state with nAv photons
    assert(nMax<101,'nMax too high for simulation');
    rho = thermalState(nMax,nAv);
    WF = real(mainWignerFromRho(rho))/64;
    qThermal = sum(WF);

    % Simulate theoretical coherent state with nAv photons
    rho = coherentState(nMax,nAv);
    WF = real(mainWignerFromRho(rho))/64;
    qCoherent = phaseAveragedDist(WF);
end

% Plot results
% edges = -30.0625:0.125:30.0625;
% xHist = -30:0.125:30;
xAxis = -20:0.125:20;

% Adjust histogram edges to discretization and plot histogram
uniq = unique(X(:));
hDisc = min(diff(uniq)); % discretization
histEdges = (min(uniq)-hDisc/2):hDisc:(max(uniq)+hDisc/2);
xHist = min(uniq):hDisc:max(uniq);
h = histogram(X,histEdges,'Normalization','probability');
h.EdgeColor = 'b';

% Plot theory curves
if nMax<101
    hold on;
    renorm = min(diff(xAxis))/hDisc;
    qThermal = qThermal/renorm;
    qCoherent = qCoherent/renorm;
    plot(xAxis,qThermal,'linewidth',2);
    plot(xAxis,qCoherent,'linewidth',2);
    legend('Measured','Thermal','Coherent');
    
    % Calculate euclidean distance from both theoretical states
    start = (length(h.Values) - length(qThermal))/2+1;
    stop = start + length(qThermal) - 1;
%     distTherm = sqrt(sum((qThermal - h.Values(start:stop)).^2));
%     distCoh = sqrt(sum((qCoherent - h.Values(start:stop)).^2));
    distTherm = [];
    distCoh = [];
else
    legend('Measured');
    [distTherm,distCoh] = deal(NaN);
end
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

