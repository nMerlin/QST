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
    Phases = mod(linspace(0,phaseperiods*360,nSamples),360)';
end

%% Step 2: Draw average photon numbers
if strcmp(type,'Thermal')
    % Draw average photon numbers from simple exponential distribution
    r = rand(nSamples,1);
    Photons = exponentialCDF(r,nPhotons,'Inverse',true);
elseif strcmp(type,'Coherent')
    % Constant average photon number
    Photons = ones(nSamples,1)*nPhotons;
end

%% Step 3: Draw quadratures
% Draw uniformly distributed probabilities to get quadratures from CDF
r = rand(nSamples,1);
X = cohCDF(r,Phases,Photons,'Inverse',true);

end

