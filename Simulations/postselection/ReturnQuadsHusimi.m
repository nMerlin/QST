function [ RetVec ] = ReturnQuadsHusimi( HusMat, MaxQuad, Resolution )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

QuadVals=-abs(MaxQuad):Resolution:abs(MaxQuad);

[XAxis,YAxis]=meshgrid(QuadVals,QuadVals);



Qx=0;
Qy=0;
Qx2=0;
Qy2=0;

HusMat=HusMat./(sum(sum(HusMat)));
Qx=sum(sum(XAxis.*HusMat));
Qy=sum(sum(YAxis.*HusMat));

FullQx2=sum(sum(XAxis.^2.*HusMat));
FullQy2=sum(sum(YAxis.^2.*HusMat));

Qx2=sum(sum((XAxis-Qx).^2.*HusMat));
Qy2=sum(sum((YAxis-Qy).^2.*HusMat));

PhotonNr=0.5*(FullQx2+FullQy2)-1; %See e.g. Carmichael 1, page 144

FullQx4=sum(sum(XAxis.^4.*HusMat));
FullQy4=sum(sum(YAxis.^4.*HusMat));


RetVec=[Qx,Qy,Qx2,Qy2,PhotonNr];

end

