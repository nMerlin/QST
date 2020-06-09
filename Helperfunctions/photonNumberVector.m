function [n] = photonNumberVector(X)
% computes a vector with photon numbers from the quadratures X, that has
% the same length as X. 
[pulses,segments,piezoSegments]=size(X);
Xres = X(:);
Xres = reshape(Xres,[pulses segments*piezoSegments]); % Consecutive values in columns
Xres = Xres - mean(mean(Xres));
ada = mean(Xres.^2)-0.5;
ada(ada<0) = 0; %set photon numbers < 0 to 0.

n = zeros(size(Xres));
for iRow = 1:pulses
    n(iRow,1:end) = ada;
end
n = n(:);
n = reshape(n,size(X));
end