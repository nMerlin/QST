function nPhotons = photonsWigner(WF, q)
%PHOTONSWIGNER Returns the mean photon number of the given Wigner function
%
%   Input Parameters:
%       WF - discretized Wigner function (normalized to 1)
%       q - values used for discretization in both dimensions
%
%   Physics:
%       <n> = 1/(4*A^2)*(<q^2>+<p^2>)-1/2
%       this formula depends on your choice of A with q=A*(aDagger + a)

A = 1/sqrt(2);

meanQ2 = mean(q.^2*sum(WF,1)');
meanP2 = mean(q.^2*sum(WF,2));
nPhotons = 1/(4*A^2)*(meanQ2 + meanP2) - 1/2;

end
