function [ QF ] = thermHusimi(q,p,nPhotons)
%THERMHUSIMI Returns Q(q,p) for a thermal state with NPHOTONS.
%   Detailed explanation goes here

disc = min(diff(q));
QF = zeros(length(q),length(p));
for iP = 1:length(p)
    QF(:,iP)=2*disc^2/pi*1/(nPhotons+1)* ...
        exp(-1/2*(q.^2+p(iP).^2)/(nPhotons+1));
end
%HF = HF./sum(sum(HF)); % Renorm (necessary due to discretization)

end

