function [nX1,nX2,nX3] = nPhotons(X1,X2,X3)
%NPHOTONS Computes the number of photons from the 3-Channel quadratures

% nX1 = var(X1(:)) - 0.5;
% nX2 = var(X2(:)) - 0.5;
% nX3 = var(X3(:)) - 0.5;

nX1 = mean(X1(:).^2) - 0.5;
nX2 = mean(X2(:).^2) - 0.5;
nX3 = mean(X3(:).^2) - 0.5;

end

