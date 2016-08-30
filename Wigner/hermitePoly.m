% HermitePoly.m by David Terr, Raytheon, 5-10-04
% https://de.mathworks.com/matlabcentral/fileexchange/4911-hermitepoly-m/content//HermitePoly.m

function hk = hermitePoly(n)
%HERMITEPOLY computes the Hermite polynomial H_n.
% Given nonnegative integer n, compute the Hermite polynomial H_n. Return
% the result as a vector whose mth element is the coefficient of x^(n+1-m).
% polyval(HermitePoly(n),x) evaluates H_n(x).

if n==0 
    hk = 1;
elseif n==1
    hk = [2 0];
else
    
    hkm2 = zeros(1,n+1);
    hkm2(n+1) = 1;
    hkm1 = zeros(1,n+1);
    hkm1(n) = 2;

    for k=2:n
        
        hk = zeros(1,n+1);

        for e=n-k+1:2:n
            hk(e) = 2*(hkm1(e+1) - (k-1)*hkm2(e));% /sqrt(2*k)
        end
        
        hk(n+1) = -2*(k-1)*hkm2(n+1);% /sqrt(2*k)
        
        if k<n
            hkm2 = hkm1;
            hkm1 = hk;
        end
        
    end
    
end