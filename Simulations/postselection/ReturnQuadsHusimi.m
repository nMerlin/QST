function [Qx,Qy,FullQx2,FullQy2,VarQx,VarQy,VarQxWigner,VarQyWigner,PhotonNr ] = ReturnQuadsHusimi( HusMat, MaxQuad, Resolution )
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

VarQxWigner = VarQx -0.5;
VarQyWigner = VarQy -0.5;
% For a gaussian signal state, then: VarHusimiSignal = VarWignerSignal +
% VarWignerCoherent, because Husimi is a convolution of Wigner and coherent
% Wigner, and the property of convolution of gaussians, see e.g. 
% https://de.wikipedia.org/wiki/Normalverteilung#Invarianz_gegen%C3%BCber_Faltung
% We assume the variance of coherent Wigner function is 0.5

PhotonNr=0.5*(FullQx2+FullQy2) -1; % See e.g. Carmichael 1, page 144. A Husimi function has variance = 1 for the vacuum state.

FullQx4=sum(sum(XAxis.^4.*HusMat));
FullQy4=sum(sum(YAxis.^4.*HusMat));




end

