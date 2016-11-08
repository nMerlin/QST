function [ PI ] = computeProjector( x, theta, maxFockState )
%COMPUTEPROJECTOR Summary of this function goes here
%   Detailed explanation goes here

PI = zeros(maxFockState+1, maxFockState+1, length(x));
for n=0:maxFockState
    PI(n+1,n+1,:) = abs(fockstate(n,x)).^2;
    for m = (n+1):maxFockState
        PI(n+1,m+1,:) = fockstate(n,x) .* fockstate(m,x) ...
            .*exp(1i*(n-m)*theta);
        PI(m+1,n+1,:) = conj(PI(n+1,m+1,:));
    end
end

end

