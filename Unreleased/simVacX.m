function [X] = simVacX(varargin)
%SIMVACX Outputs the simulated quadratures of a vacuum state

%% Validate and parse input arguments
p = inputParser;
defaultQuadsDisc = 1/2^3;
addParameter(p,'QuadsDisc',defaultQuadsDisc,@isnumeric);
defaultPhaseDisc = 1; % degrees
addParameter(p,'PhaseDisc',defaultPhaseDisc,@isnumeric);
defaultMaxQuadrature = 20;
addParameter(p,'MaxQuadrature',defaultMaxQuadrature,@isnumeric);
parse(p,varargin{:});
c = struct2cell(p.Results);
[maxquadrature,phasedisc,quadsdisc] = c{:};

%% Step 1: Discretize Wigner function and Phase Vector
[q,p] = deal(-maxquadrature:quadsdisc:maxquadrature);
WF = cohWigner(q,p,0);
PH = 0:phasedisc:(360-phasedisc);
% plotWigner(WF,'surface');

%% Step 2: Compute Integral Projections

%% Step 3: Sampling with random phase vector

end

