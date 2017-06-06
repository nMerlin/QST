function [ rho ] = thermalState( nMax, nAv )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Compute stationary state of a cavity in a thermal bath with nAv photons
a = nAv/(1+nAv); % a is exp(-beta*omega)
rho = zeros(nMax + 1);
rho(1,1) = (1-a)/(1-a^(nMax+1));
for i = 2:nMax + 1
    rho(i,i) = rho(i-1,i-1)*a;
end

end

