function [X] = simVacX(nSamples,varargin)
%SIMVACX Outputs the simulated quadratures of a vacuum state
%
%   Input Arguments:
%       nSamples: Number of quadratures to generate.

%% Validate and parse input arguments
p = inputParser;
defaultMaxQuadrature = 20;
addParameter(p,'MaxQuadrature',defaultMaxQuadrature,@isnumeric);
defaultNPhotons = 10;
addParameter(p,'NPhotons',defaultNPhotons,@isnumeric);
defaultPhaseDisc = 1; % degrees
addParameter(p,'PhaseDisc',defaultPhaseDisc,@isnumeric);
defaultQuadsDisc = 1/2^3;
addParameter(p,'QuadsDisc',defaultQuadsDisc,@isnumeric);
parse(p,varargin{:});
c = struct2cell(p.Results);
[maxquadrature,nPhotons,phasedisc,quadsdisc] = c{:};

%% Step 1: Discretize Wigner function and Phase Vector
[q,p] = deal(-maxquadrature:quadsdisc:maxquadrature);
PH = 0:phasedisc:(360-phasedisc);

%% Step 2: Cumulative distribution functions of integral projections
CDF = zeros(length(q),length(PH));
for iPhase=1:length(PH)
    WF = cohWigner(q,p,nPhotons,'Theta',PH(iPhase));
    IP = sum(WF); % Integral Projection
    CDF(:,iPhase) = cumsum(IP); % Cumulative distribution function
end

%% Step 3: Inverse transform sampling of CDF to generate Quadratures
X = zeros(nSamples,1);
Phases = unidrnd(length(PH),nSamples,1);
r = rand(nSamples,1);
for iSample=1:nSamples
    iQ = find(CDF(:,Phases(iSample))>=r(iSample),1);
    X(iSample) = q(iQ);
end

end

