function [meanPhase,meanAmp,meanPhaseBinned,meanAmpBinned,varPhase,varAmp,varPhaseBinned,varAmpBinned,meanAbsPhase,PhotonNr,PhotonNrPhase,g1,Qx,Qy,FullQx2,FullQy2,VarQx,VarQy,meanAmpQP,sigNeg] = ReturnMomentsFromP( P, sigmaP, QuadVals, varBins,filename,varargin )
%UNTITLED2 Computes the expectation values for a given
%regularized P function. 
%% Validate and parse input arguments
p = inputParser;
defaultNorm = 1;  %the standard deviation of the vacuum
addParameter(p,'Norm',defaultNorm,@isnumeric); 
defaultPlot = true;
addParameter(p,'Plot',defaultPlot,@islogical);
parse(p,varargin{:});
c = struct2cell(p.Results);
[norm,plotOption] = c{:}; 

%%
% quantum coherence?? 

%significant negativity in terms of the standard deviation 
sigVals = P./sigmaP;
sigNeg = min(min(sigVals(P < 0)));

% quadrature exp. values 
[QAxis,PAxis]=meshgrid(QuadVals,QuadVals);

P=P./(sum(sum(P)));
Qx=sum(sum(QAxis.*P));
Qy=sum(sum(PAxis.*P));

FullQx2=sum(sum(QAxis.^2.*P));%Der Wert ist zu gro? im Vergleich zum echten Wert.
FullQy2=sum(sum(PAxis.^2.*P));

VarQx=sum(sum((QAxis-Qx).^2.*P));%Der Wert ist zu gro? im Vergleich zum echten Wert.
VarQy=sum(sum((PAxis-Qy).^2.*P));%Der Wert ist zu gro? im Vergleich zum echten Wert.

PhotonNr=1/(4*norm^2)*(FullQx2+FullQy2) - 0.5;   % geht das so f?r P funktion? Einerseits muss die Normierung ber?cksichtig werden;
%Andererseits ist der Wert immer noch zu gro? im Vergleich zum echten Wert.
%

FullQx4=sum(sum(QAxis.^4.*P));
FullQy4=sum(sum(PAxis.^4.*P));

meanAmpQP=sum(sum(sqrt(QAxis.^2 + PAxis.^2).*P));

% phase and amplitude exp. values 
PhaseMatrix = atan2(PAxis,QAxis); %returns values in the closed interval [-pi,pi]
AmplMatrix = sqrt(QAxis.^2 + PAxis.^2);

% reduce the area to a circle
MaxQuad = max(QuadVals);
iCirc = find(AmplMatrix <= abs(MaxQuad));
PhaseMatrix = PhaseMatrix(iCirc);
AmplMatrix = AmplMatrix(iCirc);
P = P(iCirc);
P=P./(sum(sum(P)));
% compute expectation values 
meanPhase=sum(sum(PhaseMatrix.*P));
meanAbsPhase=sum(sum(abs(PhaseMatrix).*P));
meanAmp=sum(sum(AmplMatrix.*P));
PhotonNrPhase= 1/(4*norm^2)*sum(sum(AmplMatrix.^2.*P)) - 0.5; 
varPhase=sum(sum((PhaseMatrix-meanPhase).^2.*P));
varAmp=sum(sum((AmplMatrix-meanAmp).^2.*P));
g1 = abs(sum(sum(exp(1i*PhaseMatrix).*P)))^2;

% 


%% Sorting of P function Values and amplitude values into Phase Bins 
[N,edges,bin] = histcounts(PhaseMatrix(:),varBins);
[~,I] = sort(bin);
WFsort=P(:);
WFsort = WFsort(I);
WFOut = NaN(max(N), varBins);
for iInterval = 1 : varBins
    start = 1+sum(N(1:iInterval-1));
    stop = start+N(iInterval)-1;
    WFOut(1:N(iInterval),iInterval) = WFsort(start:stop);
end
WFPhaseBinned = sum(WFOut, 'omitnan');
WFPhaseBinned = WFPhaseBinned/sum(WFPhaseBinned, 'omitnan');
PhAxis=(0:1:varBins-1)*mean(diff(edges))+min(edges);
meanPhaseBinned = sum(PhAxis.*WFPhaseBinned, 'omitnan');
varPhaseBinned = sum(PhAxis.^2.*WFPhaseBinned, 'omitnan') - meanPhaseBinned.^2;
g1WithoutAmpBinned = abs(sum(exp(1i*PhAxis).*WFPhaseBinned))^2;

if plotOption
    fontsize2 = 15;
    bar(PhAxis,WFPhaseBinned );
    xlabel('Phase \phi = atan2(P,Q)');
    ylabel('WF(\phi)');
    text('Units','normalized','position',[0.6,0.8],'String',...
        ['<\phi> = ' num2str(meanPhase,'%.2f') char(10) 'Var(\phi) = ' ...
        num2str(varPhase,'%.2f') char(10) '<\phi>_{Binned} = ' num2str(meanPhaseBinned,'%.2f') ...
        char(10) 'Var(\phi)_{Binned} = ' ...
        num2str(varPhaseBinned,'%.2f') ],'FontSize',fontsize2);
    graphicsSettings; 
    print([filename '-Phasehistogram.png'],'-dpng');
    savefig([filename '-Phasehistogram.fig']);
end

%% Sorting of Wignerfunction Values into Amplitude Bins 
[N,edges,bin] = histcounts(AmplMatrix(:),varBins);
[~,I] = sort(bin);
WFsort=P(:);
WFsort = WFsort(I);
WFOut = NaN(max(N), varBins);
for iInterval = 1 : varBins
    start = 1+sum(N(1:iInterval-1));
    stop = start+N(iInterval)-1;
    WFOut(1:N(iInterval),iInterval) = WFsort(start:stop);
end
AmpAxis=(0:1:varBins-1)*mean(diff(edges))+min(edges);
WFAmpBinned = sum(WFOut, 'omitnan');
WFAmpBinned = WFAmpBinned/sum(WFAmpBinned, 'omitnan');
meanAmpBinned = sum(AmpAxis.*WFAmpBinned, 'omitnan');
varAmpBinned = sum(AmpAxis.^2.*WFAmpBinned, 'omitnan') - meanAmpBinned^2;

if plotOption
    bar(AmpAxis,WFAmpBinned );
    xlabel('$r = \sqrt{Q^2 + P^2}$','Interpreter','latex');
    ylabel('WF(r)');
    text('Units','normalized','position',[0.6,0.8],'String',...
        ['<r> = ' num2str(meanAmp,'%.2f') char(10) 'Var(r) = ' ...
        num2str(varAmp,'%.2f') char(10) '<r>_{Binned} = ' num2str(meanAmpBinned,'%.2f') char(10) 'Var(r)_{Binned} = ' ...
        num2str(varAmpBinned,'%.2f') ],'FontSize',fontsize2);
    graphicsSettings; 
    print([filename '-Amplitudehistogram.png'],'-dpng');
    savefig([filename '-Amplitudehistogram.fig']);
end

end