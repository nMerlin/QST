function [ WF ] = cohWigner( q, p, nPhotons )
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

sigma = 1/sqrt(2); % Select another norm if applicable
q0 = sqrt(4*sigma^2*nPhotons);
p0 = 0;

WF = zeros(length(q),length(p));
for iP = 1:length(p)
    WF(:,iP) = 1/(2*pi*sigma^2)*exp(-1/(2*sigma^2)*((q-q0).^2+ ...
        (p(iP)-p0).^2));
end
WF = WF./sum(sum(WF)); % Renorm (necessary due to discretization)

end
