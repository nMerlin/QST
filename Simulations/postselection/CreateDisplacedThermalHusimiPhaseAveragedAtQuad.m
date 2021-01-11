function [ HusFunc ] = CreateDisplacedThermalHusimiPhaseAveragedAtQuad( nThermal,r1,r2,QVal,PVal,WRes )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% thermal photon number nThermal and
% coherent radius r. It can be elliptical when chosing different r1 and r2.
% 
  
    
HusFunc = 0;
for theta = 0:0.1:2*pi; 
    Q0 = r1*cos(theta);
    P0 = r2*sin(theta);   
    HusFuncPhase=(WRes^2/pi)*(1/(nThermal+1))*exp(-0.5*((QVal-Q0).^2+(PVal-P0).^2)/(nThermal+1));
    HusFunc = HusFunc + HusFuncPhase;
end;
HusFunc=HusFunc./sum(sum(HusFunc));
    
%This defintion of the coherent state means that the quads do not correspond to
%the amplitudes, <Qx^2>=<Qy^2>=0.5 for the vacuum state, so Q
%\propto\sqrt(2)* \alpha


end

