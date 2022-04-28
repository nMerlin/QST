function [ PF ] = cohGlauberSudarshan(q,p,nPhotons,varargin)
%THERMHUSIMI Returns P(q,p) for a coherent state with NPHOTONS.
%   Formula from Schleich, Chapter 12, 12.8
%   The standard deviation, sigma=A, depends on your
%   choice of normalization q = A*(a^dagger + a). The offset is only
%   applied in the q-direction.
%
%   Input Parameters:
%       q,p - discretization vectors for the Husimi function
%       nPhotons - average number of photons

%% Validate and parse input arguments
parser = inputParser;
defaultP0 = sqrt(2*nPhotons); % Change that for another normalization!
addParameter(parser,'P0',defaultP0,@isnumeric);
defaultQ0 = 0;
addParameter(parser,'Q0',defaultQ0,@isnumeric);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[p0,q0] = c{:};

%%

p = p';
PF_q = dirac(q-q0);
idx = PF_q == Inf;
PF_q(idx) = 1; 
PF_p = dirac(p-p0);
idx = PF_p == Inf;
PF_p(idx) = 1; 
PF = PF_q + PF_p;



end