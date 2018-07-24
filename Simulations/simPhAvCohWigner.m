function WF = simPhAvCohWigner(q,p,nPhotons,nTheta)
%SIMPHAVWIGNER Returns a phase-averaged coherent Wigner function
%
%   q,p: Discretized quadrature axes
%   nPhotons: Average number of photons of the coherent state
%   nTheta: Number of phases to use for the averaging.

WF = zeros(length(q),length(p));
for i = 0:(nTheta-1)
    WF = WF + cohWigner(q,p,nPhotons,'Theta',i*360/nTheta);
end
WF = WF./sum(sum(WF));

end

