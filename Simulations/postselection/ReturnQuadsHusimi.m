function [Qx,Qy,FullQx2,FullQy2,VarQx,VarQy,PhotonNr ] = ReturnQuadsHusimi( HusMat, MaxQuad, Resolution )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

QuadVals=-abs(MaxQuad):Resolution:abs(MaxQuad);

[XAxis,YAxis]=meshgrid(QuadVals,QuadVals);




HusMat=HusMat./(sum(sum(HusMat)));
Qx=sum(sum(XAxis.*HusMat));
Qy=sum(sum(YAxis.*HusMat));

FullQx2=sum(sum(XAxis.^2.*HusMat));
FullQy2=sum(sum(YAxis.^2.*HusMat));

VarQx=sum(sum((XAxis-Qx).^2.*HusMat));
VarQy=sum(sum((YAxis-Qy).^2.*HusMat));

PhotonNr=0.5*(FullQx2+FullQy2-1); %PhotonNr=0.5*(FullQx2+FullQy2)-1; See e.g. Carmichael 1, page 144??

FullQx4=sum(sum(XAxis.^4.*HusMat));
FullQy4=sum(sum(YAxis.^4.*HusMat));




end

