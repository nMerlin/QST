function PSHusimi = PostselectedHusimiDisplacedThermalVac(Q2,P2,nThermal,Q0,P0,Transmission)
%Calculates Wigner function at position (Q1,P1) for mixing displaced thermal light of photon number n and
%vacuum, then performing postselection on the Q outcome in channel 2
%
% Transmission gives transmission coefficient of the beam splitter, lokk in
% incerse BS trafo
%
% 
% 
% 

P2res=0.1;
P2min=-20;
P2max=20;


QuadVals=P2min:P2res:P2max;

[xx,yy] = meshgrid(QuadVals,QuadVals);

[ModQ1,ModP1,ModQ2,ModP2]=InverseBSTrafo( xx,yy,Q2,P2,Transmission );


    %%%%DISPLACED THERMAL AND VACUUM
    PSHusimi=(1/P2res^2) *(CreateDisplacedThermalHusimiAtQuad(nThermal,Q0,P0,ModQ1,ModP1,P2res).*CreateCoherentHusimiAtQuad( 0,0,ModQ2,ModP2,P2res ));
    %%%%DISPLACED THERMAL AND VACUUM
    