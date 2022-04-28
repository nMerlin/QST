function [ HusFunc ] = CreateDisplacedThermalHusimiAtQuad( n,Q0,P0,QVal,PVal,WRes )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    HusFunc=(WRes^2/pi)*(1/(n+1))*exp(-0.5*((QVal-Q0).^2+(PVal-P0).^2)/(n+1));
    
%This defintion of the coherent state means that the quads do not correspond to
%the amplitudes, <Qx^2>=<Qy^2>=0.5 for the vacuum state, so Q
%\propto\sqrt(2)* \alpha


end

