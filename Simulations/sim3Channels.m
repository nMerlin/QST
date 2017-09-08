function [X1,X2,X3] = sim3Channels()
%SIMCROSSCORR Simulates a 3-Channel thermal state measurement.

%% Constants
% Original Number of Signal Photons
meanN = 30;

% Beamsplitter coefficients
t1 = sqrt(1/30); % transmitted beam goes to channel 3
r1 = sqrt(1-t1^2); % reflected beam gets split again for channels 1 and 2
t2 = 1/sqrt(2);
r2 = sqrt(1-t2^2);

% Amplitude of vacuum noise
noiseampl = 1/sqrt(2);

% Modulation frequencies in Hz
freq1 = 50;
freq2 = 25;
freq3 = 0.5;

% Buffers for trigger
nPeriodsPerSeg = 1.6;
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
% Gaussian Noise to simulate vacuum in each homodyne channel
Ns1 = reshape(randn(N,1)*noiseampl,[nPulses,nPieces,nSegments]);
Ns2 = reshape(randn(N,1)*noiseampl,[nPulses,nPieces,nSegments]);
Ns3 = reshape(randn(N,1)*noiseampl,[nPulses,nPieces,nSegments]);

% Gaussian Noise to simulate vacuum mixed at beamsplitters for copying
avac1 = reshape(randn(N,1)*noiseampl,[nPulses,nPieces,nSegments]);
avac2 = reshape(randn(N,1)*noiseampl,[nPulses,nPieces,nSegments]);

% Amplitude Noise according to Bose-Einstein distribution (thermal state)
thermN = reshape(inverseBoseEinstein(rand(N,1)/(meanN+1),meanN), ...
    [nPulses,nPieces,nSegments]);

% Compute Quadratures from Noise and Phase
X1 = Ns1 + r1*t2*sqrt(thermN) .* sin(phase1) - t1*t2*avac1 + r2*avac2;
X2 = Ns2 + r1*r2*sqrt(thermN) .* sin(phase2) - t1*r2*avac1 - t2*avac2;
X3 = Ns3 + t1*sqrt(thermN) .* sin(phase3) + r1*avac1;

end

function N = inverseBoseEinstein(y,meanN)
    N = (log(y)-log(1/(meanN+1)))/log(meanN/(meanN+1));
end

