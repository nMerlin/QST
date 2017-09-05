function [cc,X1,X2] = simCrossCorr()
%SIMCROSSCORR Simulates the cross-correlation between two homodyne
%detection channels.

% Constants
nPulses = 1000;
nPieces = 400;
nSegments = 1;
N = nPulses*nPieces*nSegments;
theta = linspace(pi/2,1.6*2*pi,N)';
ampl = 5;
noiseampl = 1/sqrt(2);

% Uniformly distributed phase values
phi = rand(N,1)*2*pi;

% Gaussian Noise to simulate vacuum
Ns1 = randn(N,1)*noiseampl;
Ns2 = randn(N,1)*noiseampl;

% Measured Quadratures
X1 = Ns1 + ampl*sin(phi);
X2 = Ns2 + ampl*sin(phi+theta);
X1 = reshape(X1,[nPulses,nPieces,nSegments]);
X2 = reshape(X2,[nPulses,nPieces,nSegments]);

% Cross-Correlation
cc = X1.*X2;

end

