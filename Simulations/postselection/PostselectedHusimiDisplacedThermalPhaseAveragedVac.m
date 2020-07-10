function PSHusimi = PostselectedHusimiDisplacedThermalPhaseAveragedVac(Q2,P2,n,r1,r2,Transmission,MaxQuad,Resolution)
%Calculates Husimi function at position (Q1,P1) for mixing vacuum and phase 
%averaged displaced thermal light of total photon number n and
%radius r (or different radii r1 and r2), then performing postselection on the Q outcome in channel 2 (Q2,
%P2 correspond to the postselection channel); they are scaled with sqrt(2)
%because the postselection arm is divided into two channels. 
%
% Transmission gives transmission coefficient of the beam splitter, look in
% inverse BS trafo
%
% 
% 
% 

QuadVals=-abs(MaxQuad):Resolution:abs(MaxQuad);

[xx,yy] = meshgrid(QuadVals,QuadVals);

[ModQ1,ModP1,ModQ2,ModP2]=InverseBSTrafo( xx,yy,sqrt(2)*Q2,sqrt(2)*P2,Transmission );

% n is the total photon number, consisting of nThermal and nCoherent, which is given by radius r.
% get nThermal:
nThermal = n - 0.5*mean([r1^2 r2^2]);

%%%%DISPLACED THERMAL AND VACUUM
H_PDT = CreateDisplacedThermalHusimiPhaseAveragedAtQuad(nThermal,r1,r2,ModQ1,ModP1,Resolution);
H_PDT = H_PDT / sum(sum(H_PDT));
Hvac = CreateCoherentHusimiAtQuad( 0,0,ModQ2,ModP2,Resolution );
Hvac = Hvac / sum(sum(Hvac));

PSHusimi=(1/Resolution^2) *H_PDT.*Hvac;
%%%%DISPLACED THERMAL AND VACUUM
    
end