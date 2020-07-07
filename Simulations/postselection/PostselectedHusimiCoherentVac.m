function PSHusimi = PostselectedHusimiCoherentVac(Q2,P2,n,Transmission)
%Calculates Wigner function at position (Q1,P1) for mixing coherent light of photon number n and
%vacuum, then performing postselection on the Q outcome in channel 2
% Using n is a stub. Giving the detailed quadratures in both directions
% would be better.
% Transmission gives transmission coefficient of the beam splitter
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


    %%%%COHERENT AND VACUUM
    PSHusimi=(1/P2res^2) *(CreateCoherentHusimiAtQuad(sqrt(2*n),0,ModQ1,ModP1,P2res).*CreateCoherentHusimiAtQuad( 0,0,ModQ2,ModP2,P2res ));
    %%%%COHERENT AND VACUUM
    