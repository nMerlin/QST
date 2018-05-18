function minVar = compute3ChLimit(nX1,nX2,nX3)
%COMPUTE3CHLIMIT Theoretical minimum variance of 3Ch measurement.
%
% Note: nX3 is the photon number of the target channel.

n = nX1 + nX2 + nX3;
nt = nX3;
nps = nX1 + nX2;
minVar = (1+n+nt)/(2*(1+nps));

end

