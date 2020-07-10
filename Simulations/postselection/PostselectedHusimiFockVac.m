function PSHusimi = PostselectedHusimiFockVac(Q2,P2,n,Transmission)
%Calculates Wigner function at position (Q1,P1) for mixing a Fock state of photon number n and
%vacuum, then performing postselection on the Q outcome in channel 2 (Q2,
%P2 correspond to the postselection channel); they are scaled with sqrt(2)
%because the postselection arm is divided into two channels. 
%
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

[ModQ1,ModP1,ModQ2,ModP2]=InverseBSTrafo( xx,yy,sqrt(2)*Q2,sqrt(2)*P2,Transmission );


    %%%%FOCK AND VACUUM
    PSHusimi=(1/P2res^2) *(CreateFockHusimiAtQuad(n,ModQ1,ModP1,P2res).*CreateCoherentHusimiAtQuad( 0,0,ModQ2,ModP2,P2res ));
    %%%%FOCK AND VACUUM
    