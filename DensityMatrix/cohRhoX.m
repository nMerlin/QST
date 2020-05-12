function [ rhoX ] = cohRhoX( nAv,varargin )
%COHERENTSTATE Coherent state density matrix in x,y basis
% for psi_x : Davidovic 1993, doi = {10.1088/0305-4470/26/19/037}, eq. 2
%   Parameters
%   NAV - Average photon number of the coherent state

%% Validate and parse input arguments
parser = inputParser;
defaultXY = -20:0.125:20;
addParameter(parser,'XY',defaultXY,@isvector);
defaultNorm = 1/sqrt(2); % sqrt([q,p]/(2*i))
addParameter(parser,'Norm',defaultNorm,@isnumeric);
defaultP0 = sqrt(2*nAv); % Change that for another normalization!
addParameter(parser,'P0',defaultP0,@isnumeric);
defaultQ0 = 0;
addParameter(parser,'Q0',defaultQ0,@isnumeric);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[norm,p0,q0,xy] = c{:};

%% prepare x, y axis
[x, y] = deal(xy);
nX = length(x);

%% make p and q matrix
xMatrix = zeros(nX,nX);
for index = 1:nX
    xMatrix(index,:) = x;
end
yMatrix = xMatrix';

%% make rho
psiX = psi(xMatrix,q0,p0);
psiY = psi(yMatrix,q0,p0);
rhoX = psiX.*conj(psiY);


%% Wavefunction Psi_x = <x|alpha>
    function [psiX] = psi(x,q0,p0)
        psiX = pi^0.25 * exp( -0.5*(x - q0).^2 + 1i*p0*x);
    end

end