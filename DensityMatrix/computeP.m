function [ P ] = computeP( maxFockStates )
%COMPUTEP quadrature operator P in fock basis

P = zeros(maxFockStates+1, maxFockStates+1);
for m = 1:maxFockStates+1
    for n = 1:maxFockStates+1
        P(m,n) = sqrt(n/2) * (m==n-1) - sqrt((n+1)/2) * (m==n+1);
    end
end

end

