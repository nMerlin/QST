function [n] = probThermPhotonsInverse(nProb,nAverage)
%PROBTHERMPHOTONS Returns probability of measuring n thermal photons.
n = (log(nProb)+log(nAverage+1))/(log(nAverage)-log(nAverage+1));
end

