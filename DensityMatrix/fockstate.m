function fval = fockstate(n,x)
%HERMITEGAUSSIAN computes the number state wavefunction <n,x>.

%% Simple method with hermite Polynomials
% This method is numerically less stable because of factorial(n)
%fval = polyval(hermitePoly(n),x).*(pi^(-0.25)).* ...
%    (1/sqrt((2^n)*factorial(n))).*exp(-0.5*(x.^2));

%% Recursion
% See, for example, https://doi.org/10.1088/1361-6404/aa9584 with alpha = 1
% Takes very long for high n
% if n == 0
%     fval = pi^(-0.25)*exp(-0.5*x.^2);
% elseif n == 1
%     fval = pi^(-0.25)*exp(-0.5*x.^2).*x*sqrt(2);
% else
%     fval = sqrt(2/n)*x.*fockstate(n-1,x)-sqrt((n-1)/n)*fockstate(n-2,x);
% end

%% Iteration
% Iterative derivation using the recursion formula
F_2 = pi^(-0.25)*exp(-0.5*x.^2);
F_1 = pi^(-0.25)*exp(-0.5*x.^2).*x*sqrt(2);
F_0 = 0;
if n == 0
    fval = F_2;
elseif n == 1
    fval = F_1;
else
    for iN = 2:n
        F_0 = sqrt(2/iN)*x.*F_1-sqrt((iN-1)/iN)*F_2;
        F_2 = F_1;
        F_1 = F_0;
    end
    fval = F_0;
end

end
