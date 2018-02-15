function [ WF ] = cohWigner( q, p, nPhotons, varargin )
%COHWIGNER Returns W(q,p) for a coherent state with NPHOTONS photons
%
%   The Wigner function of a coherent state is a Gaussian function with an
%   offset. However, the standard deviation, sigma=A, depends on your
%   choice of normalization q = A*(a^dagger + a). The offset is only
%   applied in the q-direction.
%
%   Input Parameters:
%       q,p - discretization vectors for the Wigner function
%       nPhotons - average number of photons

%% Validate and parse input arguments
parser = inputParser;
defaultTheta = 0; % Rotation angle in degrees
addParameter(parser,'Theta',defaultTheta,@isnumeric);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[theta] = c{:};

sigma = 1/sqrt(2); % Select another norm if applicable
q0 = 0;
p0 = sqrt(2*nPhotons); % Change that for another normalization!

WF = zeros(length(q),length(p));
for iP = 1:length(p)
    qrot = q*cosd(theta)-p(iP)*sind(theta);
    prot = q*sind(theta)+p(iP)*cosd(theta);
    WF(:,iP) = 1/(2*pi*sigma^2)*exp(-1/(2*sigma^2)* ...
        ((qrot-q0).^2+((prot-p0)).^2));
end
WF = WF./sum(sum(WF)); % Renorm (necessary due to discretization)

end
