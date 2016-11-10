function [ rho ] = coherentState( nMax, Ns )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Coherent state density matrix in the Fock (number) state representation
rho =zeros(nMax+1,1);
% In this program, all matrices are diagonal
if (Ns <100)
    for I=1:(nMax+1)
        rho(I) =exp(-Ns) *Ns^(I-1) /factorial(I-1);
    end
else% Use Stirling's approximation
    for I=1:(nMax+1)
        rho(I) =exp(I-1-Ns) *(Ns/(I-1))^(I-1) /sqrt(2*pi*(I-1));
    end
end

rho = diag(rho);

end

