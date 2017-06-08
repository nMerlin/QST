function [ nAv ] = plotPhAvData( X )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Remove offset
X = X(:);
X = X - mean(mean(X));

% Mean photon number
nAv = 0.5*(2*mean(X.^2))-0.5;

% Simulate theoretical thermal state with nAv photons
nMax = max(ceil(3*nAv),30);
rho = thermalState(nMax,nAv);
WF = real(mainWignerFromRho(rho))/64;
qThermal = sum(WF);

% Simulate theoretical coherent state with nAv photons
rho = coherentState(nMax,nAv);
WF = real(mainWignerFromRho(rho))/64;
qCoherent = phaseAveragedDist(WF);

% Plot results
edges = -20.0625:0.125:20.0625;
histogram(X,edges,'Normalization','probability');
hold on;
plot(-20:0.125:20,qThermal,'linewidth',2);
plot(-20:0.125:20,qCoherent,'linewidth',2);

end

