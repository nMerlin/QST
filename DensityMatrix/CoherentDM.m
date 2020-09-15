function [ DM ] = CoherentDM( X,P,MaxN )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% X,P: Offset in phase space
DM=zeros(MaxN+1,MaxN+1);

Alpha=(X+1i*P)/sqrt(2);

for n1=0:MaxN
    for n2=0:MaxN
        DM(n1+1,n2+1)=exp(-abs(Alpha).^2) .*Alpha.^n1 .*conj(Alpha).^n2/sqrt(factorial(n1)*factorial(n2));
    end
end

end 