function HG = hermiteGaussian(n,x)
%HERMITEGAUSSIAN computes the number state wavefunction <x|n>.
HG=polyval(hermitePoly(n),x).*(pi^(-0.25)).*(1/sqrt((2^n)*factorial(n))).*exp(-0.5*(x.^2));
end
