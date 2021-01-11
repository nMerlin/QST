function [ HF ] = fockHusimi(q,p,nPhotons,varargin)
%THERMHUSIMI Returns Q(q,p) for a coherent state with NPHOTONS.
%   Formula from Schleich, Chapter 12, 12.9
%   The standard deviation, sigma=A, depends on your
%   choice of normalization q = A*(a^dagger + a). The offset is only
%   applied in the q-direction.
%
%   Input Parameters:
%       q,p - discretization vectors for the Husimi function
%       nPhotons - average number of photons, must be integer

%% Validate and parse input arguments
parser = inputParser;
defaultNorm = 1/sqrt(2); % sqrt([q,p]/(2*i))
addParameter(parser,'Norm',defaultNorm,@isnumeric);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[norm] = c{:};

%%
disc = mean(diff(q));
HF = zeros(length(q),length(p));
for iP = 1:length(p)
    absAlphaSquared = (q.^2 + p(iP).^2)/(2*norm)^2;
    HF(:,iP)=disc^2 *1/pi* absAlphaSquared.^(nPhotons) / factorial(nPhotons) .* exp(-absAlphaSquared);
end
HF = HF./sum(sum(HF)); % Usually not necessary

end