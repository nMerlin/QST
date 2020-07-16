function [Qx,Qy,FullQx2,FullQy2,VarQx,VarQy,PhotonNr ] = ReturnQuadsWigner( WF, MaxQuad, Resolution )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

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




end

