function [X] = simMultiChannels(nSamples,varargin)
%SIMVACX Outputs the simulated quadratures of a vacuum state
%
%   Input Arguments:
%       nSamples: Number of quadratures to generate.

%% Validate and parse input arguments
p = inputParser;
defaultLoadFile = '';
addParameter(p,'LoadFile',defaultLoadFile,@isstr);
defaultMaxQuadrature = 20;
addParameter(p,'MaxQuadrature',defaultMaxQuadrature,@isnumeric);
defaultNPhotons = 10;
addParameter(p,'NPhotons',defaultNPhotons,@isnumeric);
defaultPhaseDisc = 1; % degrees
addParameter(p,'PhaseDisc',defaultPhaseDisc,@isnumeric);
defaultPhasePeriods = 1;
addParameter(p,'PhasePeriods',defaultPhasePeriods,@isnumeric);
defaultPhotonDisc = 1;
addParameter(p,'PhotonDisc',defaultPhotonDisc,@isnumeric);
defaultQuadsDisc = 1/2^3;
addParameter(p,'QuadsDisc',defaultQuadsDisc,@isnumeric);
defaultSaveFile = '';
addParameter(p,'SaveFile',defaultSaveFile,@isstr);
defaultType = 'Thermal';
addParameter(p,'Type',defaultType,@isstr);
parse(p,varargin{:});
c = struct2cell(p.Results);
[loadfile,maxquadrature,nPhotons,phasedisc,phaseperiods,photondisc, ...
    quadsdisc,savefile,type] = c{:};

if strcmp(loadfile,'')
    %% Step 1: Discretize Wigner function, Phase Vector and Photon Numbers
    [q,p] = deal(-maxquadrature:quadsdisc:maxquadrature);
    PH = 0:phasedisc:(360-phasedisc); % Discretized phases
    PV = 0:photondisc:5*nPhotons; % Discretized photons (photon vector)
    nPH = length(PH);
    nPV = length(PV);

    %% Step 2: Cumulative distribution functions of integral projections
    CDF = zeros(length(q),length(PH),length(PV));
    parfor iPhase=1:nPH
        for iPhotons=1:nPV
            WF = cohWigner(q,p,PV(iPhotons),'Theta',PH(iPhase));
            IP = sum(WF); % Integral Projection
            CDF(:,iPhase,iPhotons) = cumsum(IP); % Cumulative distribution func
        end
        iPhase
    end

    % Save CDF matrix and other parameters
    if not(strcmp(savefile,''))
        save(savefile,'q','PH','PV','nPH','nPV','CDF');
    end
else
    load(loadfile);
end

%% Step 3: Inverse transform sampling of CDF to generate Quadratures
X = zeros(nSamples,1);

%% Step 3.1: Draw phases
if strcmp(type,'Thermal')
    % Draw uniformly distributed phases
    iPhases = unidrnd(nPH,nSamples,1);
elseif strcmp(type,'Coherent')
    Phases = mod(linspace(0,phaseperiods*360,nSamples),360);
    [~,iPhases] = histc(Phases, ...
        [-Inf (phasedisc/2):phasedisc:(359-phasedisc/2) Inf]);
end

%% Step 3.2: Draw photon numbers
probs = probThermPhotons(PV,nPhotons);
probs = probs./sum(probs);
photonCDF = cumsum(probs);
nPhotonCDF = length(photonCDF);
if strcmp(type,'Thermal')
    % Draw Bose-Einstein distributed photon numbers and find their
    % corresponding indices in PV
    [~,iPhotons] = histc(rand(nSamples,1), ...
        [-Inf interp1(1:nPhotonCDF,photonCDF,0.5+(1:nPhotonCDF-1)) Inf]);
elseif strcmp(type,'Coherent')
    [~,iPhotons] = histc(ones(nSamples,1)*nPhotons, ...
        [-Inf interp1(1:nPV,PV,0.5+(1:nPV-1)) Inf]);
end

%% Step 3.3: Draw quadratures
% Draw uniformly distributed probabilities to get quadratures from CDF
r = rand(nSamples,1);
for iSample=1:nSamples
    iQ = find(CDF(:,iPhases(iSample),iPhotons(iSample))>=r(iSample),1);
    X(iSample) = q(iQ);
%     X(iSample) = interp1(CDF(:,Phases(iSample),iPhotons(iSample)), ...
%         q,r(iSample));
end

end

