function [thetaEven,thetaOdd] = alternatingPiezoSegmentsTheta(theta)
% takes a phase vector theta that is devided in piezo segments and creates a new
% vector containing only the even piezo segments, and another vector
% containing the odd piezosegments. 
[seg,piezoSeg] = size(theta);
thetaEven = zeros([seg floor(piezoSeg/2)]);
thetaOdd = zeros([seg ceil(piezoSeg/2)]);
for j = 1:piezoSeg
   if rem(j,2) == 0       
       thetaEven(:,j/2) = theta(:,j); %select only even piezosegments
   else
       thetaOdd(:,(j+1)/2) = theta(:,j); %select only odd piezosegments
   end    
end
end