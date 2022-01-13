function [Xnew] = flipPiezoSegments(X)
% takes a vector X that is devided in piezo segments and flips every other
% piezo segment in time.
[pulses,seg,piezoSeg] = size(X);
Xnew = reshape(X,[pulses*seg piezoSeg]);
for j = 1:piezoSeg
   if rem(j,2) == 1       
        Xnew(:,j) = flip(Xnew(:,j));
   end    
end
Xnew = reshape(Xnew,[pulses seg piezoSeg]);
end

