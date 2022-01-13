function [Xeven,Xodd] = alternatingPiezoSegments(X)
% takes a vector X that is devided in piezo segments and creates a new
% vector containing only the even piezo segments, and another vector
% containing the odd piezosegments. 
[pulses,seg,piezoSeg] = size(X);
X = reshape(X,[pulses*seg piezoSeg]);
Xeven = zeros([pulses*seg floor(piezoSeg/2)]);
Xodd = zeros([pulses*seg ceil(piezoSeg/2)]);
for j = 1:piezoSeg
   if rem(j,2) == 0       
       Xeven(:,j/2) = X(:,j); %select only even piezosegments
   else
       Xodd(:,(j+1)/2) = X(:,j); %select only odd piezosegments
   end    
end
Xeven = reshape(Xeven,[pulses seg floor(piezoSeg/2)]);
Xodd = reshape(Xodd,[pulses seg ceil(piezoSeg/2)]);
end