function [nProb] = probThermPhotons(n,nAverage)
%PROBTHERMPHOTONS Returns probability of measuring n thermal photons.
nProb = 1/(nAverage+1)*(nAverage/(nAverage+1)).^n;
end

