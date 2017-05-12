function [ Q ] = computeQ( maxFockStates )
%COMPUTEQ quadrature operator Q in fock basis

Q = zeros(maxFockStates+1, maxFockStates+1);
for m = 1:maxFockStates+1
    for n = 1:maxFockStates+1
        Q(m,n) = sqrt(n/2) * (m==n-1) + sqrt((n+1)/2) * (m==n+1);
    end
end

end

