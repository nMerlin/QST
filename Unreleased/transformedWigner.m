function [QF,Q,Q2,P,P2] = transformedWigner(q,p)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

QF = zeros(length(q),length(p));
for iP = 1:length(p)
    QF(:,iP) = 1/pi*exp(-1/2*(q.^2+p(iP).^2));
end
QF = QF./sum(sum(QF)); % Renorm (necessary due to discretization)

%% Expectation values
Q = sum(QF*q);
P = sum(QF'*p);
Q2 = sum(QF*(q-Q).^2);
P2 = sum(QF'*(p-P).^2);

end

