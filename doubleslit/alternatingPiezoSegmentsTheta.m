function [thetaEven,thetaOdd,thetaKorr] = alternatingPiezoSegmentsTheta(theta, piezoSign)
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

%change timing sequenz in each piezosegment for the reverse driving piezo
if piezoSign == -1
    thetaOdd = flip(thetaOdd);  
else
    thetaEven = flip(thetaEven);
end

thetaKorr = zeros(size(theta));

for k = 1: piezoSeg
    if rem(k,2) == 0       
      thetaKorr(:,k)= thetaEven(:,k/2);  
   else
       thetaKorr(:,k)= thetaOdd(:,(k+1)/2); 
   end    
end
thetaKorr = reshape(thetaKorr,[123752,513]);
end