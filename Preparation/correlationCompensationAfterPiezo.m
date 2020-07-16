function [Xrem] = correlationCompensationAfterPiezo(X)
%removes the correlations of quadratures that are already piezoshaped. 
   [pulses,stuff,segments] = size(X);
   Xrem = reshape(X,[pulses,stuff*segments]);
   Xrem = correlationCompensation(Xrem);
   Xrem = Xrem(1:end-1,:);  
   Xrem = reshape(Xrem,[pulses-1,stuff,segments]);
end