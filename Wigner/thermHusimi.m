function [ HF ] = thermHusimi(q,p,nPhotons)
%THERMHUSIMI Returns Q(q,p) for a thermal state with NPHOTONS.
%   Detailed explanation goes here

HF = zeros(length(q),length(p));
for iP = 1:length(p)
    HF(:,iP) = 1/pi*1/(nPhotons+1)*exp(-(q.^2+p(iP).^2)/(nPhotons+1));
end
HF = HF./sum(sum(HF)); % Renorm (necessary due to discretization)

end

