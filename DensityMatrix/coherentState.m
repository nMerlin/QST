function [ rho ] = coherentState( nMax, nAv )
%COHERENTSTATE Coherent state density matrix in Fock basis
%   Paramters:
%   NMAX - Size of the density matrix minus 1
%   NAV - Average photon number of the coherent state

% Coherent state density matrix in the Fock (number) state representation
alpha = sqrt(nAv);
rho = zeros(nMax+1,1);
for n = 1:(nMax+1)
    for m = 1:(nMax+1)
        rho(n,m) = alpha^(n+m)/sqrt(factorial(n)*factorial(m))*...
            exp(-abs(alpha)^2);
    end
end

end
