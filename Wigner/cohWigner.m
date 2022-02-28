function [ WF ] = cohWigner( q, p, nPhotons, varargin )
%COHWIGNER Returns W(q,p) for a coherent state with NPHOTONS photons
%
%   The Wigner function of a coherent state is a Gaussian function with an
%   offset. However, the standard deviation, norm=A, depends on your
%   choice of normalization q = A*(a^dagger + a). The offset is only
%   applied in the q-direction.
%
%   Input Parameters:
%       q,p - discretization vectors for the Wigner function
%       nPhotons - average number of photons

%% Validate and parse input arguments
parser = inputParser;
defaultNorm = 1/sqrt(2); % sqrt([q,p]/(2*i))
addParameter(parser,'Norm',defaultNorm,@isnumeric);
defaultP0 = sqrt(2*nPhotons); % Change that for another normalization!
addParameter(parser,'P0',defaultP0,@isnumeric);
defaultQ0 = 0;
addParameter(parser,'Q0',defaultQ0,@isnumeric);
defaultTheta = 0; % Rotation angle in degrees
addParameter(parser,'Theta',defaultTheta,@isnumeric);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[norm,p0,q0,theta] = c{:};

WF = zeros(length(q),length(p));
for iP = 1:length(p)
    qrot = q*cosd(theta)-p(iP)*sind(theta);
    prot = q*sind(theta)+p(iP)*cosd(theta);
    WF(:,iP) = 1/(2*pi*norm^2)*exp(-1/(2*norm^2)* ...
        ((qrot-q0).^2+((prot-p0)).^2));
end
WF = WF./sum(sum(WF)); % Renorm (necessary due to discretization)

end
