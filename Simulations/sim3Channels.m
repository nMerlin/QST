function [X1,X2,X3] = sim3Channels()
%SIMCROSSCORR Simulates a 3-Channel thermal state measurement.

%% Constants
% Number of Photons
nPhot1 = 10;
nPhot2 = 10;
nPhot3 = 10;
noiseampl = 1/sqrt(2);

% Modulation frequencies in Hz
freq1 = 50;
freq2 = 25;
freq3 = 0.5;

% Buffers for trigger
nPeriodsPerSeg = 2.0;
upBuf = 0.2;
lowBuf = 0.25;

% Sizing of arrays
nPulses = 999;
nPieces = 406;
nSegments = 157;
N = nPulses*nPieces*nSegments;
nEdgePieces = round(nPieces/(1-upBuf-lowBuf));
lenEdges = nPulses*nEdgePieces;

%% Simulation of Piezo Modulation and Triggers
t = (1:lenEdges*nSegments)*pi/lenEdges;
t = reshape(t,[nPulses,nEdgePieces,nSegments]);
startLow = round(nEdgePieces*lowBuf);
startUp = round(nEdgePieces*upBuf);
tLow = t(:,startLow:startLow+nPieces-1,1:2:nSegments);
tUp = t(:,startUp:startUp+nPieces-1,2:2:nSegments);
t = zeros(nPulses,nPieces,nSegments);
itup = 1; itlow = 1;
for k = 1:size(tLow,3)+size(tUp,3)
    if mod(k,2)==1
        t(:,:,k) = tLow(:,:,itlow);
        itlow = itlow+1;
    else
        t(:,:,k) = tUp(:,:,itup);
        itup = itup+1;
    end
end
phaseScaling = nPeriodsPerSeg*2*pi/(1-lowBuf-upBuf);
phase1 = abs(sawtooth(t))*phaseScaling;
phase2 = abs(sawtooth(t*freq2/freq1))*phaseScaling;
phase3 = abs(sawtooth(t*freq3/freq1))*phaseScaling;

%% Add random phase noise
% Uniformly distributed phase values
randPhase = reshape(rand(N,1)*2*pi,[nPulses,nPieces,nSegments]);
phase1 = phase1 + randPhase;
phase2 = phase2 + randPhase;
phase3 = phase3 + randPhase;

%% Compute Quadratures
% Gaussian Noise to simulate vacuum
Ns1 = reshape(randn(N,1)*noiseampl,[nPulses,nPieces,nSegments]);
Ns2 = reshape(randn(N,1)*noiseampl,[nPulses,nPieces,nSegments]);
Ns3 = reshape(randn(N,1)*noiseampl,[nPulses,nPieces,nSegments]);

% Compute Quadratures from Noise and Phase
X1 = Ns1 + sqrt(nPhot1) * sin(phase1);
X2 = Ns2 + sqrt(nPhot2) * sin(phase2);
X3 = Ns3 + sqrt(nPhot3) * sin(phase3);

end

