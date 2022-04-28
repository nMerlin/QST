function [meanPhase,meanAmp,meanPhaseBinned,meanAmpBinned,varPhase,varAmp,varPhaseBinned,varAmpBinned,meanAbsPhase,PhotonNr,g1WithAmp,g1WithoutAmp,g1WithoutAmpBinned,g1WithAmpNorm] = ReturnPhaseAndAmplitudeWigner( WF, MaxQuad, Resolution, varBins,filename,varargin )
%UNTITLED2 Computes the expectation values of Phase and Amplitude, derived from quadratures, for a given
%Wigner function. The maximum quadrature value MaxQuad and their Resolution
%need to match those used for computing the Wigner function. 
%% Validate and parse input arguments
p = inputParser;
defaultPlot = true;
addParameter(p,'Plot',defaultPlot,@islogical);
parse(p,varargin{:});
c = struct2cell(p.Results);
[plotOption] = c{:};

%%
QuadVals=-abs(MaxQuad):Resolution:abs(MaxQuad);

[QAxis,PAxis]=meshgrid(QuadVals,QuadVals);
PhaseMatrix = atan2(PAxis,QAxis); %returns values in the closed interval [-pi,pi]
AmplMatrix = sqrt(QAxis.^2 + PAxis.^2);

% reduce the area to a circle
iCirc = find(AmplMatrix <= abs(MaxQuad));
PhaseMatrix = PhaseMatrix(iCirc);
AmplMatrix = AmplMatrix(iCirc);
WF=WF./(sum(sum(WF)));
WF = WF(iCirc);
WF=WF./(sum(sum(WF)));
% compute expectation values 
meanPhase=sum(sum(PhaseMatrix.*WF));
meanAbsPhase=sum(sum(abs(PhaseMatrix).*WF));
meanAmp=sum(sum(AmplMatrix.*WF));
PhotonNr=0.5*sum(sum(AmplMatrix.^2.*WF)) - 0.5; 
varPhase=sum(sum((PhaseMatrix-meanPhase).^2.*WF));
varAmp=sum(sum((AmplMatrix-meanAmp).^2.*WF));
g1WithAmp = abs(sum(sum(AmplMatrix.*exp(1i*PhaseMatrix).*WF)))^2;
g1WithAmpNorm = g1WithAmp/abs(sum(sum(AmplMatrix.*WF)))^2; %normalizing
g1WithoutAmp = abs(sum(sum(exp(1i*PhaseMatrix).*WF)))^2;

%% Sorting of Wignerfunction Values and amplitude values into Phase Bins 
[N,edges,bin] = histcounts(PhaseMatrix(:),varBins);
[~,I] = sort(bin);
WFsort=WF(:);
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
%g1WithAmpBinned = abs(sum(AmplMatrix.*exp(1i*PhAxis).*WFPhaseBinned))^2;
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
WFsort=WF(:);
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