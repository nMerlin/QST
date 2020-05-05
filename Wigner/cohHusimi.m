function [ HF ] = cohHusimi(q,p,nPhotons,varargin)
%THERMHUSIMI Returns Q(q,p) for a coherent state with NPHOTONS.
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
defaultNorm = 1/sqrt(2); % sqrt([q,p]/(2*i))
addParameter(parser,'Norm',defaultNorm,@isnumeric);
defaultP0 = sqrt(2*nPhotons); % Change that for another normalization!
addParameter(parser,'P0',defaultP0,@isnumeric);
defaultQ0 = 0;
addParameter(parser,'Q0',defaultQ0,@isnumeric);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[norm,p0,q0] = c{:};

%%
disc = mean(diff(q));
HF = zeros(length(q),length(p));
for iP = 1:length(p)
    HF(:,iP)=disc^2 *1/pi* exp(-((q - q0).^2 + (p(iP) - p0).^2)/(2*norm)^2);
end
HF = HF./sum(sum(HF)); % Usually not necessary

end