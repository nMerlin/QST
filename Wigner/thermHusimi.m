function [ QF ] = thermHusimi(q,p,nPhotons,disc)
%THERMHUSIMI Returns Q(q,p) for a thermal state with NPHOTONS.
%   Detailed explanation goes here

QF = zeros(length(q),length(p));
for iP = 1:length(p)
    QF(:,iP)=disc^2*1/(pi*(nPhotons+1))*exp(-(q.^2+p(iP).^2)/(nPhotons+1));
end
%HF = HF./sum(sum(HF)); % Usually not necessary

end

