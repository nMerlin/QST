function [A_c, Var_q,Var_p, n_c] = evaluatePostselection(nTotal,r1,r2,MaxQuad,Resolution,varargin)
% makes simulations of postselection outcomes.
% Right now it only works for phase averaged displaced thermal states.
% nTotal: total photon number of input signal state before all beam
% splitters
% r1,r2 coherent radius of signal state before all beam
% splitters
% MaxQuad, Resolution: maximum quadrature value and quadrature Resolution
% Transmission: Intensity ratio of target photon number / total photon
% number
% A_ps: postselection radius after beam splitters
% Thickness: range for postselection values around A_ps
% PhaseAverage: True if the postselection quadratures should assume all
% phases

%% Validate and parse input arguments
p = inputParser;
defaultTransmission = 0:0.1:1; %chose 0.23 for present experiment
addParameter(p,'Transmission',defaultTransmission,@isvector);
defaultA_ps = 0:1:15;
addParameter(p,'A_ps',defaultA_ps,@isvector);
defaultThickness = 0;
addParameter(p,'Thickness',defaultThickness,@isnumeric);
defaultPhaseAverage = false;
addParameter(p,'PhaseAverage',defaultPhaseAverage,@islogical);
defaultPhaseUncertainty = false;
addParameter(p,'PhaseUncertainty',defaultPhaseUncertainty,@islogical);
defaultNDisc = 100;
addParameter(p,'NDisc',defaultNDisc,@isnumeric);
defaultState = 'PhaseAveragedDisplacedThermal';
addParameter(p,'State',defaultState,@isstr);
parse(p,varargin{:});
c = struct2cell(p.Results);
[A_ps,nDisc,phaseAverage,phaseUncertainty,state,thickness,Transmission] = c{:};

if phaseAverage
    theta = 0:0.1:2*pi;
else    
    theta = 0;
end

if phaseUncertainty
    lowerLim = -pi/2;
    upperLim = pi/2;
    Res = (upperLim-lowerLim)/10;    
    thetaUnc = lowerLim:Res:upperLim;
else
    thetaUnc = 0;
end

[A_c, Var_q, Var_p, n_c] = deal(zeros(length(Transmission),length(A_ps)));

for i = 1:length(Transmission)
    for j = 1: length(A_ps)
        A_psVector = A_ps(j)-thickness/2:0.1:A_ps(j)+thickness/2;       
        for k = 1:length(theta) 
            PSHusimi = 0;
            for m = 1:length(A_psVector)  
                for u = 1:length(thetaUnc)
                    switch state
                        case 'PhaseAveragedDisplacedThermal'                    
                            PSHusimi = PSHusimi + PostselectedHusimiDisplacedThermalPhaseAveragedVac(A_psVector(m)*cos(theta(k)+thetaUnc(u)),...
                                A_psVector(m)*sin(theta(k)+thetaUnc(u)),nTotal,r1,r2,Transmission(i),MaxQuad,Resolution);
                        case 'Thermal'
                            PSHusimi = PSHusimi + PostselectedHusimiThermalVac(A_psVector(m)*cos(theta(k)+thetaUnc(u)),...
                                A_psVector(m)*sin(theta(k)+thetaUnc(u)),nTotal,Transmission(i),MaxQuad,Resolution);
                        case 'DisplacedThermal'
      
                            PSHusimi = PSHusimi + PostselectedHusimiDisplacedThermalVac(A_psVector(m)*cos(theta(k)+thetaUnc(u)),...
                                A_psVector(m)*sin(theta(k)+thetaUnc(u)),nTotal,5,0,Transmission(i),MaxQuad,Resolution);
                    end
                end
            end
                PSHusimi = PSHusimi / sum(sum(PSHusimi));
%                 [WF] = WignerFromHusimi(PSHusimi,'PQ',-abs(MaxQuad):Resolution:abs(MaxQuad));
%                 [Qx,Qy,~,~,VarQxWigner,VarQyWigner,PhotonNr ] = ReturnQuadsWigner( WF, MaxQuad, Resolution );
                [Qx,Qy,~,~,~,~,VarQxWigner,VarQyWigner,PhotonNr] = ReturnQuadsHusimi( PSHusimi, MaxQuad, Resolution );
                A_c(i,j) = A_c(i,j) + sqrt(Qx^2+Qy^2)/length(theta);
                Var_q(i,j) = Var_q(i,j) + VarQxWigner /length(theta);
                Var_p(i,j) = Var_p(i,j) + VarQyWigner /length(theta);
                n_c(i,j) = n_c(i,j) + PhotonNr/length(theta);               
        end
    end
end

%% Plot Amplitude
figure(1);
if length(Transmission) == 1
    plot(A_ps,A_c)
    ylabel('A_c')
    title(['Transmission = ' num2str(Transmission)]);
else
    surf(A_ps,Transmission,A_c);
    shading('flat');
    ylabel('Transmission');
    zlabel('A_c');
end
xlabel('A_{ps}');
graphicsSettings;
savefig(['A_c-nTotal-' num2str(nTotal) '-r1-' num2str(r1) '-r2-' num2str(r2) '-maxQuad-' ...
    num2str(MaxQuad) '-Res-' num2str(Resolution) '-phaseUncertainty-' num2str(phaseUncertainty) ...
    '-nDisc-' num2str(nDisc) '-state-' num2str(state) '.fig']);

%% Plot Variance
figure(2);
if length(Transmission) == 1
    plot(A_ps,Var_q,A_ps,Var_p);
    legend('Var_{Q}','Var_{P}','location','best');
    ylabel('Variance (Wigner)')
    xlabel('A_{ps}');
    %ylim([0.5 5]);
    title(['Transmission = ' num2str(Transmission)]);
    graphicsSettings;
else
    subplot(1,2,1);
    surf(A_ps,Transmission,Var_q);
    shading('flat');
    ylabel('Transmission');
    zlabel('Variance_Q (Wigner)');
    xlabel('A_{ps}');
    %zlim([0.5 1]);
    graphicsSettings;
    
    subplot(1,2,2);
    surf(A_ps,Transmission,Var_p);
    shading('flat');
    ylabel('Transmission');
    xlabel('A_{ps}');
    zlabel('Variance_P (Wigner)');
    %zlim([0.5 1]);
    graphicsSettings;
end

savefig(['Var-nTotal-' num2str(nTotal) '-r1-' num2str(r1) '-r2-' num2str(r2) '-maxQuad-' ...
    num2str(MaxQuad) '-Res-' num2str(Resolution) '-phaseUncertainty-' num2str(phaseUncertainty) ...
    '-nDisc-' num2str(nDisc) '-state-' num2str(state) '.fig']);

%% Plot Photon Number
figure(3);
if length(Transmission) == 1
    plot(A_ps,n_c)
    ylabel('postselected Photon Number')
    title(['Transmission = ' num2str(Transmission)]);
else
    surf(A_ps,Transmission,n_c);
    shading('flat');
    ylabel('Transmission');
    zlabel('postselected Photon Number');
end
xlabel('A_{ps}');
graphicsSettings;
savefig(['n_c-nTotal-' num2str(nTotal) '-r1-' num2str(r1) '-r2-' num2str(r2) '-maxQuad-' ...
    num2str(MaxQuad) '-Res-' num2str(Resolution) '-phaseUncertainty-' num2str(phaseUncertainty) ...
    '-nDisc-' num2str(nDisc) '-state-' num2str(state) '.fig']);





end