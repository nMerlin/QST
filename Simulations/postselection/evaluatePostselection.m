function [A_c, Var_c, n_c] = evaluatePostselection(nTotal,r1,r2,MaxQuad,Resolution,varargin)
% makes simulations of postselection outcomes.
% Right now it only works for phase averaged displaced thermal states.
% nTotal: total photon number of input signal state
% r1,r2 coherent radius of signal state 
% MaxQuad, Resolution: maximum quadrature value and quadrature Resolution
% Transmission: Intensity ratio of target photon number / total photon
% number
% A_ps: postselection radius 
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
parse(p,varargin{:});
c = struct2cell(p.Results);
[phaseAverage,thickness,Transmission] = c{:};

if phaseAverage
    theta = 0:0.1:2*pi;
else    
    theta = 0;
end

[A_c, Var_c, n_c] = deal(zeros(length(Transmission),length(A_ps)));

for i = 1:length(Transmission)
    for j = 1: length(A_ps)
        A_psVector = A_ps(j)-thickness/2:0.1:A_ps(j)+thickness/2;       
        for k = 1:length(theta) 
            PSHusimi = 0;
            for m = 1:length(A_psVector)                     
                PSHusimi = PSHusimi + PostselectedHusimiDisplacedThermalPhaseAveragedVac(A_psVector(m)*cos(theta(k)),...
                    A_psVector(m)*sin(theta(k)),nTotal,r1,r2,Transmission(i),MaxQuad,Resolution);                
            end
                PSHusimi = PSHusimi / sum(sum(PSHusimi));
                [Qx,Qy,~,~,VarQx,VarQy,PhotonNr] = ReturnQuadsHusimi( PSHusimi, MaxQuad, Resolution );
                A_c(i,j) = A_c(i,j) + sqrt(Qx^2+Qy^2)/length(theta);
                Var_c(i,j) = Var_c(i,j) + mean([VarQx VarQy])/length(theta);
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
    num2str(MaxQuad) '-Res-' num2str(Resolution) '.fig']);

%% Plot Variance
figure(2);
if length(Transmission) == 1
    plot(A_ps,Var_c)
    ylabel('Variance')
    title(['Transmission = ' num2str(Transmission)]);
else
    surf(A_ps,Transmission,Var_c);
    shading('flat');
    ylabel('Transmission');
    zlabel('Variance');
end
xlabel('A_{ps}');
graphicsSettings;
savefig(['Var_c-nTotal-' num2str(nTotal) '-r1-' num2str(r1) '-r2-' num2str(r2) '-maxQuad-' ...
    num2str(MaxQuad) '-Res-' num2str(Resolution) '.fig']);

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
    num2str(MaxQuad) '-Res-' num2str(Resolution) '.fig']);





end