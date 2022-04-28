function [ HusFunc ] = CreateFockHusimiAtQuad( n,QVal,PVal,HRes )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


HusFunc=(HRes^2/pi).* exp(-0.5*(QVal.^2+PVal.^2)).*((0.5*(QVal.^2+PVal.^2)).^n)./factorial(n);
    
%This defintion of the coherent state means that the quads do not correspond to
%the amplitudes, <Q^2>=0.5


end

