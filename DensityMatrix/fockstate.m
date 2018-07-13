function fockstate = fockstate(n,x)
%HERMITEGAUSSIAN computes the number state wavefunction <n,x>.
fockstate=polyval(hermitePoly(n),x).*(pi^(-0.25)).* ...
    (1/sqrt((2^n)*factorial(n))).*exp(-0.5*(x.^2));
end
