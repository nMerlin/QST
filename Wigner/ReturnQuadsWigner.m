function [Qx,Qy,FullQx2,FullQy2,VarQx,VarQy,PhotonNr,meanAmpQP ] = ReturnQuadsWigner( WF, MaxQuad, Resolution )
%UNTITLED2 Computes the expectation values of quadratures for a given
%Wigner function. The maximum quadrature value MaxQuad and their Resolution
%need to match those used for computing the Wigner function. 

QuadVals=-abs(MaxQuad):Resolution:abs(MaxQuad);

[XAxis,YAxis]=meshgrid(QuadVals,QuadVals);

WF=WF./(sum(sum(WF)));
Qx=sum(sum(XAxis.*WF));
Qy=sum(sum(YAxis.*WF));

FullQx2=sum(sum(XAxis.^2.*WF));
FullQy2=sum(sum(YAxis.^2.*WF));

VarQx=sum(sum((XAxis-Qx).^2.*WF));
VarQy=sum(sum((YAxis-Qy).^2.*WF));

PhotonNr=0.5*(FullQx2+FullQy2) - 0.5; 

FullQx4=sum(sum(XAxis.^4.*WF));
FullQy4=sum(sum(YAxis.^4.*WF));

meanAmpQP=sum(sum(sqrt(XAxis.^2 + YAxis.^2).*WF));


end

