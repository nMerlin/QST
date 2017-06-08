function [ nAv ] = plotPhAvData( X, filename )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Remove offset
X = X(:);
X = X - mean(mean(X));

% Mean photon number
nAv = mean(X.^2)-0.5;

% Simulations only if nMax<101
nMax = max(ceil(3*nAv),30);
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
edges = -20.0625:0.125:20.0625;
xAxis = -20:0.125:20;
h = histogram(X,edges,'Normalization','probability');
if nMax<101
    hold on;
    plot(xAxis,qThermal,'linewidth',2);
    plot(xAxis,qCoherent,'linewidth',2);
    legend('Measured','Thermal','Coherent');
else
    legend('Measured');
end
xlabel('Q');
ylabel('Probability');
title(['Phase-averaged quantum state (n=',num2str(nAv),')']);

% Adjust axis limits
minX = -2*fwhm(xAxis,h.Values);
xlim([minX -minX]);
hold off;

% Save figure
print(strcat('PHAV-',filename,'-',num2str(nAv),'photons', '-dpng');

end

