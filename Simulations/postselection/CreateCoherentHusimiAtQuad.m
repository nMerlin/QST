function [ HusFunc ] = CreateCoherentHusimiAtQuad( X0,P0,XVal,PVal,Hres )
%UNTITLED2 Summary of this function goes here
%   AlphaR0: Coherent Offset (RealAxis), AlphaI0: Coherent Offset (Imag
%   Axis)

    HusFunc=(Hres^2/(pi))*exp(-0.5*((XVal-X0).^2+(PVal-P0).^2));
    
%


end

