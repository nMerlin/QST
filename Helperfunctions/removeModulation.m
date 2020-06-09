function [X1rem,X2rem,X3rem,n1,n2,n3] = removeModulation(X1,X2,X3)
%rescales the Xi from a 3 channel quadrature measurement according to their
%time dependent photon numbers in order to remove photon number
%fluctuations.

n1 = photonNumberVector(X1);
n2 = photonNumberVector(X2);
n3 = photonNumberVector(X3);

% only devide by n when it is > 0.1 to avoid deviding by too small values
X1rem = X1;
X2rem = X2;
X3rem = X3;
X1rem(n1>1) = X1(n1>1).*sqrt(mean(n1(:))./n1(n1>1));
X2rem(n2>1) = X2(n2>1).*sqrt(mean(n2(:))./n2(n2>1));
X3rem(n3>1) = X3(n3>1).*sqrt(mean(n3(:))./n3(n3>1));



end