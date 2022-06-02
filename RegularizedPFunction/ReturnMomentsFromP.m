function [meanPhase,meanAmp,meanPhaseBinned,meanAmpBinned,varPhase,circularVariance1,circularVariance2,circVar1Err,circVar2Err,varAmp,varPhaseBinned,varAmpBinned,PhotonNr,PhotonNrBinned,g1,sigNeg,phaseErr,ampErr,varPhaseErr,varAmpErr,PhotonNrErr,maxQ] = ReturnMomentsFromP( P, sigmaP, QuadVals, varBins,filename,varargin )
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
% make sure P is normalized  
sumP = sum(sum(P));
P=P./sumP;
sigmaP=sigmaP./sumP;

%significant negativity in terms of the standard deviation 
sigVals = P./sigmaP;
sigNeg = min(min(sigVals(P < 0)));

% quadrature exp. values 
[QAxis,PAxis]=meshgrid(QuadVals,QuadVals);

% phase and amplitude exp. values 
PhaseMatrix = atan2(PAxis,QAxis); %returns values in the closed interval [-pi,pi]
AmplMatrix = sqrt(QAxis.^2 + PAxis.^2);

% get the q0 value were P has its maximum
[~,maxI] = max(P(floor(length(P)/2),:));
maxQ = QuadVals(maxI);

% reduce the area to a circle
MaxQuad = max(QuadVals);
iCirc = find(AmplMatrix <= abs(MaxQuad));
PhaseMatrix = PhaseMatrix(iCirc);
AmplMatrix = AmplMatrix(iCirc);
P = P(iCirc);
sumP = sum(sum(P));
P=P./(sum(sum(P)));
sigmaP = sigmaP(iCirc)/sumP;

% compute expectation values 
meanPhase=sum(sum(PhaseMatrix.*P));
meanAmp=sum(sum(AmplMatrix.*P));
varPhase=sum(sum((PhaseMatrix-meanPhase).^2.*P));
varAmp=sum(sum((AmplMatrix-meanAmp).^2.*P));
phaseErr = sqrt(sum(sigmaP.^2.*abs(PhaseMatrix).^2));
phaseSquareErr = sqrt(sum(sigmaP.^2.*abs(PhaseMatrix).^4));
varPhaseErr = sqrt( phaseSquareErr^2 + (-2*meanPhase)^2 * phaseErr^2);
ampErr = sqrt(sum(sigmaP.^2.*abs(AmplMatrix).^2));
ampSquareErr = sqrt(sum(sigmaP.^2.*abs(AmplMatrix).^4));
varAmpErr = sqrt( ampSquareErr^2 + (-2*meanAmp)^2 * ampErr^2);
g1 = abs(sum(sum(exp(1i*PhaseMatrix).*P)))^2;
PhotonNr= 1/(4*norm^2)*sum(sum(AmplMatrix.^2.*P)) - 0.5; %this seems to be too high?%
PhotonNrErr = 1/(4*norm^2)*sqrt(sum(sum(AmplMatrix.^4.*sigmaP.^2)));
circularVariance1 = 1 - sqrt( sum(sum(cos(PhaseMatrix).*P))^2 + sum(sum(sin(PhaseMatrix).*P))^2  );
circularVariance2 = 1 -  sum(sum(cos(PhaseMatrix).*P))^2 + sum(sum(sin(PhaseMatrix).*P))^2  ;
circVar2Err = sqrt( 4*sum(sum(cos(PhaseMatrix).*P))^2 * sum(sum(abs(cos(PhaseMatrix)).^2.*sigmaP.^2)) +...
    4*sum(sum(sin(PhaseMatrix).*P))^2 * sum(sum(abs(sin(PhaseMatrix)).^2.*sigmaP.^2)));
circVar1Err = 1/2 *1/(sqrt( sum(sum(cos(PhaseMatrix).*P))^2 + sum(sum(sin(PhaseMatrix).*P))^2  )) * circVar2Err;


%% Compute marginal P(phase), by summarizing all points where the phase is in a bin
varBins = varBins + mod(varBins,2);
PhAxis=(-pi:2*pi/varBins:pi); %makes sure the middle bin lies around zero. 
d = mean(diff(PhAxis));
P_ph = zeros(size(PhAxis));
for j = 1:length(PhAxis)
    P_ph(j) = sum(P(PhaseMatrix > PhAxis(j) - d/2 & PhaseMatrix < PhAxis(j) + d/2 ));
end
P_ph = P_ph/sum(P_ph, 'omitnan');
%P_ph(P_ph<0)=0;

meanPhaseBinned = sum(sum(PhAxis.*P_ph, 'omitnan'));
varPhaseBinned = sum(sum(PhAxis.^2.*P_ph, 'omitnan')) - meanPhaseBinned.^2;

if plotOption
    fontsize2 = 15;
    bar(PhAxis,P_ph );
    xlabel('Phase \phi = atan2(P,Q)');
    ylabel('P(\phi)'); 
    text('Units','normalized','position',[0.6,0.8],'String',...
        ['<\phi> = ' num2str(meanPhase,'%.2f') char(10) 'Var(\phi) = ' ...
        num2str(varPhase,'%.2f') char(10) '<\phi>_{Binned} = ' num2str(meanPhaseBinned,'%.2f') ...
        char(10) 'Var(\phi)_{Binned} = ' ...
        num2str(varPhaseBinned,'%.2f') ],'FontSize',fontsize2);
    graphicsSettings; 
    print([filename '-Phasehistogram.png'],'-dpng');
    savefig([filename '-Phasehistogram.fig']);
    close();
end

%% Compute marginal P(r)
[~,edges,~] = histcounts(AmplMatrix(:),varBins); %choses good bins for the distribution 
AmpAxis=(0:1:varBins-1)*mean(diff(edges))+min(edges);
P_r = zeros(size(AmpAxis));
for j = 1:length(AmpAxis)
    P_r(j) = sum(P(AmplMatrix > edges(j) & AmplMatrix < edges(j+1)));
end
P_r = P_r/sum(P_r, 'omitnan');
%P_r(P_r<0)=0;

meanAmpBinned = sum(AmpAxis.*P_r, 'omitnan');
varAmpBinned = sum(AmpAxis.^2.*P_r, 'omitnan') - meanAmpBinned^2;
PhotonNrBinned= 1/(4*norm^2)*sum(AmpAxis.^2.*P_r) - 0.5;

if plotOption
    bar(AmpAxis,P_r );
    xlabel('$r = \sqrt{Q^2 + P^2}$','Interpreter','latex');
    ylabel('P(r)');
    text('Units','normalized','position',[0.6,0.8],'String',...
        ['<r> = ' num2str(meanAmp,'%.2f') char(10) 'Var(r) = ' ...
        num2str(varAmp,'%.2f') char(10) '<r>_{Binned} = ' num2str(meanAmpBinned,'%.2f') char(10) 'Var(r)_{Binned} = ' ...
        num2str(varAmpBinned,'%.2f') ],'FontSize',fontsize2);
    graphicsSettings; 
    print([filename '-Amplitudehistogram.png'],'-dpng');
    savefig([filename '-Amplitudehistogram.fig']);
    close();
end

end