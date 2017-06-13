function WF = thermWigner( q, p, nPhotons )
%THERMWIGNER Returns W(q,p) for a thermal state with NPHOTONS photons
%
%   Important: Normalization A = 1/sqrt(2) for quadrature operator q =
%   A*(a^dagger + a) is used.

WF = zeros(length(q),length(p));
for iP = 1:length(p)
    WF(:,iP) = 1/pi*1/(2*nPhotons+1)*exp(-(q.^2+p(iP).^2)/(2*nPhotons+1));
end
WF = WF./sum(sum(WF)); % Renorm (necessary due to discretization)

end

