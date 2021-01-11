function [ Q1In,P1In,Q2In,P2In ] = InverseBSTrafo( Q1Out,P1Out,Q2Out,P2Out,IntTransmission )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Perform inverse beamsplitter transformation
%P1Out,P2Out,Q1Out,Q2Out are the quadrature at the BS output ports
%ModP1In,ModP2In,ModQ1In,ModQ2In are the quadratures at the input ports
%IntTransmission is the intensity transmission ration of the beamsplitter

T=sqrt(IntTransmission);
R=sqrt(1-IntTransmission);

Q1In=T*Q1Out+R*Q2Out;
P1In=T*P1Out+R*P2Out;
Q2In=R*Q1Out-T*Q2Out;
P2In=R*P1Out-T*P2Out;

end

