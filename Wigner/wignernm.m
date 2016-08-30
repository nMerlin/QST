function Wnm = wignernm(n,m,qval,pval,qintstep,minq,maxq)
%WIGNERNM evaluates the Wigner function at (qval,pval)
qrange=minq:qintstep:maxq;
Wnm=(1/(2*pi))*sum(qintstep*(hermiteGaussian(n,qval+(0.5*qrange)).*hermiteGaussian(m,qval-(0.5*qrange)).*exp(-1i*pval.*qrange)));
end