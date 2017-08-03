
function [q2Exp, q2squareExp, Var] = WignerTest2Ch(nTarget, nPS, QPS, region)

%This function calculates expectation values for 2 channels.
%One is the postselection chaneel with quadartures q1, p1. 
%The target channel has q2, p2.

%Norm ist the factor in the relation between the quadratures and the ladder
%operators: q = Norm*(a^{+} + a), p = Norm*i*(a^{+} - a)
%typical values are 1/sqrt(2) or 1/2.
Norm = 1/sqrt(2);

[q1, p1, q2, p2] = deal(linspace(-10,10,100));
%nTarget = 7.5;
%nPS = 7.5;
n = nTarget + nPS;
r = sqrt(nPS / n);
t = sqrt(nTarget / n);  

WF = Wigner2( q1, p1, q2, p2,n ,r, t , Norm);

width = 0.5; %selection width
 

switch region
    case 'fullcircle'
        R = zeros(length(q1),length(p1));
        for iQ1 = 1:length(q1)
            R(iQ1,:) = sqrt(q1(iQ1)^2+p1.^2);
        end

        [iQ1Select,iP1Select] = find(R>QPS-width & R<QPS+width);

        WF2 = sum(WF( iQ1Select, iP1Select,: ,:),1); % integrate over q1
        WF2 = sum(WF2( :, iP1Select ,:,:),2);% integrate over p1
        
        
    case 'Qline'
        %Integrate over P1:
        WF2 = sum(WF,2);

        %Select on Q1 = Qps +- width
        I = find(q1 <= QPS + width & q1 >= QPS - width);
        WF2 = sum(WF2( I, :, : ,:),1); 
end

%calculate expecatationvalues
WF2 = reshape(WF2, [100 100]);
WF2 = sum(WF2,2); %integrate over p32
WF2 = WF2/sum(WF2); %for some reason, renormalization is necessary.
q2Exp = q2*WF2;
q2squareExp = (q2.*q2)*WF2;
Var = q2squareExp - q2Exp;
end


function WF = Wigner2( q1, p1, q2, p2,nPhotons,r, t ,Norm)
%THERMWIGNER Returns W(q,p) for a thermal state with NPHOTONS photons
%
%   Important: Normalization A = 1/sqrt(2) for quadrature operator q =
%   A*(a^dagger + a) is used.


WF = zeros(length(q1),length(p1),length(q2),length(p2));
for iQ2 = 1:length(q2)
    for iP2 = 1:length(p2)
        for iP1 = 1:length(p1)
            WF(:,iP1,iQ2,iP2) = 1/pi^2 *1/(2*nPhotons+1)*...
                exp((-(r*p2(iP2)-t*p1(iP1)).^2 - (r*q2(iQ2) - t*q1).^2)/(2*Norm)^2).*...
                exp((-(r*p1(iP1)+t*p2(iP2)).^2 - (t*q2(iQ2) + r*q1).^2)/((2*nPhotons+1)*(2*Norm)^2));
        end
    end
end

WF = WF./sum(sum(sum(sum(WF)))); % Renorm (necessary due to discretization)

end

