function X = simQST(nSamples,varargin)
%SIMVACX Outputs the simulated quadratures of a vacuum state
%
%   Input Arguments:
%       nSamples: Number of quadratures to generate.

%% Validate and parse input arguments
p = inputParser;
defaultNPhotons = 10;
addParameter(p,'NPhotons',defaultNPhotons,@isnumeric);
defaultPhasePeriods = 1;
addParameter(p,'PhasePeriods',defaultPhasePeriods,@isnumeric);
defaultType = 'Thermal';
addParameter(p,'Type',defaultType,@isstr);
parse(p,varargin{:});
c = struct2cell(p.Results);
[nPhotons,phaseperiods,type] = c{:};

%% Step 1: Draw phases
if strcmp(type,'Thermal')
    % Draw uniformly distributed phases
    Phases = rand(nSamples,1)*360;
elseif strcmp(type,'Coherent')
    % Pre-defined phases (e.g. by piezo modulation)
    Phases = mod(linspace(0,phaseperiods*360,nSamples),360);
end

%% Step 2: Draw photon numbers
PV = 0:0.1:5*nPhotons;
probs = probThermPhotons(PV,nPhotons);
probs = probs./sum(probs);
photonCDF = cumsum(probs);
nPhotonCDF = length(photonCDF);
if strcmp(type,'Thermal')
    % Draw Bose-Einstein distributed photon numbers and find their
    % corresponding indices in PV
    [~,iPhotons] = histc(rand(nSamples,1), ...
        [-Inf interp1(1:nPhotonCDF,photonCDF,0.5+(1:nPhotonCDF-1)) Inf]);
    Photons = PV(iPhotons);
elseif strcmp(type,'Coherent')
    [~,iPhotons] = histc(ones(nSamples,1)*nPhotons, ...
        [-Inf interp1(1:nPV,PV,0.5+(1:nPV-1)) Inf]);
    Photons = PV(iPhotons);
end

%% Step 3: Draw quadratures
% Draw uniformly distributed probabilities to get quadratures from CDF
r = rand(nSamples,1);
X = cohCDF(r,Phases,Photons','Inverse',true);

end

