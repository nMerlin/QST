function [QF,Q,P,Q2,P2] = postselectedHusimi(Qps,Pps,n,Transmission)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Constants
resolution = 0.1;
q = -20:resolution:20;
p = -20:resolution:20;

%% Postselection on Husimi-Q
QF = zeros(length(q),length(p));
T=sqrt(Transmission);
R=sqrt(1-Transmission);
for iQ=1:length(q)
   for iP=1:length(p)
        %Perform inverse beamsplitter transformation
        %P1,P2,Q1,Q2 are the quadrature at the BS output ports
        %ModP1,ModP2,ModQ1,ModQ2 are the quadratures at the input ports
        Q1 = Qps;
        P1 = Pps;
        Q2 = q(iQ);
        P2 = p(iP);
        
        % Beta
        ModQ1=R*Q1+T*Q2;
        ModP1=R*P1+T*P2;
        
        % Gamma
        ModQ2=R*Q2-T*Q1;
        ModP2=R*P2-T*P1;
        
        QF(iQ,iP) = thermHusimi(ModQ1,ModP1,n,resolution)* ...
            thermHusimi(ModQ2,ModP2,0,resolution);
   end
end
QF = QF./sum(sum(QF));

%% Expectation values
q = q'; p = p';
Q = sum(QF*q);
P = sum(QF'*p);
Q2 = sum(QF*(q-Q).^2);
P2 = sum(QF'*(p-P).^2);

end

