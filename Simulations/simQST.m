function [X1,X2,X3,fastPhases] = simQST(nPulses,nPieces,nSegments,varargin)
%SIMQST Outputs the simulated quadratures of a quantum state tomography
%experiment
%
%   Input Arguments:
%       nSamples: Number of quadratures to generate.

%% Validate and parse input arguments
p = inputParser;
defaultChannels = 1;
addParameter(p,'Channels',defaultChannels,@isnumeric);
defaultFrequencies = [0.5 50 0];
addParameter(p,'Frequencies',defaultFrequencies,@isvector);
defaultNPhotons = 10;
addParameter(p,'NPhotons',defaultNPhotons,@isnumeric);
defaultPhasePeriods = 1;
addParameter(p,'PhasePeriods',defaultPhasePeriods,@isnumeric);
defaultR1 = 1/sqrt(3);
addParameter(p,'R1',defaultR1,@isnumeric);
defaultR2 = 1/sqrt(2);
addParameter(p,'R2',defaultR2,@isnumeric);
defaultType = 'Thermal';
addParameter(p,'Type',defaultType,@isstr);
parse(p,varargin{:});
c = struct2cell(p.Results);
[channels,frequencies,nPhotons,phaseperiods,r1,r2,type] = c{:};

%% Step 1: Draw phases (source)
nSamples = nPulses * nPieces;
if strcmp(type,'Thermal')
    % Draw uniformly distributed phases
    Phases1 = rand(nSamples,nSegments)*360;
elseif strcmp(type,'Coherent')
    % Pre-defined phases (e.g. by piezo modulation)
    Phases1 = mod(linspace(0,phaseperiods*360,nSamples),360)';
    Phases1 = repmat(Phases1,1,nSegments);
end

%% Step 2: Draw average photon numbers (source)
if strcmp(type,'Thermal')
    % Draw average photon numbers from simple exponential distribution
    r = rand(nSamples,nSegments);
    Photons1 = exponentialCDF(r,nPhotons,'Inverse',true);
elseif strcmp(type,'Coherent')
    % Constant average photon number
    Photons1 = ones(nSamples,nSegments)*nPhotons;
end

%% Step 3: Model beamsplitter effects
switch channels
    case 2
        Photons2 = Photons1 * (1-r2^2);
        Photons1 = Photons1 * r2^2;
        deltaPhases = repmat(linspace(0,phaseperiods*360,nSamples)', ...
            1,nSegments);
        Phases2 = Phases1 + deltaPhases;
    case 3
        Photons3 = Photons1 * (1-r1^2) * (1-r2^2);
        Photons2 = Photons1 * (1-r1^2) * r2^2;
        Photons1 = Photons1 * r1^2;
        % Index channels according to sorted frequencies
        [~,ind] = sort(frequencies);
        factor1 = frequencies(ind(end))/frequencies(ind(end-1));
        factor2 = frequencies(ind(end))/frequencies(ind(end-2));
        % Compute phase shifts
        Phases = zeros(nSamples,nSegments,3);
        Phases(:,:,3) = repmat(linspace(0,phaseperiods*360,nSamples)', ...
            1,nSegments); % fastest modulation
        Phases(:,:,2) = reshape(linspace(0, ...
            360*phaseperiods/factor1*nSegments,nSamples*nSegments), ...
            [nSamples,nSegments]); % second fastest
        Phases(:,:,1) = reshape(linspace(0, ...
            360*phaseperiods/factor2*nSegments,nSamples*nSegments), ...
            [nSamples,nSegments]); % third fastest
        % Add phaseshifts to corresponding channels
        Phases3 = Phases1 + Phases(:,:,ind==3);
        Phases2 = Phases1 + Phases(:,:,ind==2);
        Phases1 = Phases1 + Phases(:,:,ind==1);
end

%% Step 4: Draw quadratures
% Draw uniformly distributed probabilities to get quadratures from CDF
[X2,X3] = deal([]);
X1 = cohCDF(rand(nSamples,nSegments),Phases1, ...
    Photons1,'Inverse',true);
X1 = reshape(X1,[nPulses,nPieces,nSegments]);
if channels > 1
    X2 = cohCDF(rand(nSamples,nSegments),Phases2, ...
        Photons2,'Inverse',true);
    X2 = reshape(X2,[nPulses,nPieces,nSegments]);
end
if channels > 2
    X3 = cohCDF(rand(nSamples,nSegments),Phases3, ...
        Photons3,'Inverse',true);
    X3 = reshape(X3,[nPulses,nPieces,nSegments]);
end

end

