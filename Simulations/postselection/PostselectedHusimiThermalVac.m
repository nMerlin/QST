function PSHusimi = PostselectedHusimiThermalVac(Q2,P2,n,Transmission,MaxQuad,Resolution)
%Calculates Wigner function at position (Q1,P1) for mixing thermal light of photon number n and
%vacuum, then performing postselection on the Q outcome in channel 2.(Q2,
%P2 correspond to the postselection channel); they are scaled with sqrt(2)
%because the postselection arm is divided into two channels. 
%
% Transmission gives transmission coefficient of the beam splitter
%
% 
% 
% 


QuadVals=-abs(MaxQuad):Resolution:abs(MaxQuad);

[xx,yy] = meshgrid(QuadVals,QuadVals);

[ModQ1,ModP1,ModQ2,ModP2]=InverseBSTrafo( xx,yy,sqrt(2)*Q2,sqrt(2)*P2,Transmission );


    %%%%THERMAL AND VACUUM
    PSHusimi=(1/Resolution^2) *(CreateThermalHusimiAtQuad(n,ModQ1,ModP1,Resolution).*CreateCoherentHusimiAtQuad( 0,0,ModQ2,ModP2,Resolution ));
    %%%%THERMAL AND VACUUM
    