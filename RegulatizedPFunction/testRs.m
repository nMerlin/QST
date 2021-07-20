%testing P function for different Rs
%Normalization: If the standard deviation of the vacuum is set to 1 
selX2 = sqrt(2)*selX;
maxQuad = 20;
Res = 0.25;
maxX = 20;
XStep = 0.25;
phiStep = 0.1;
dispstat('','init','timestamp','keepthis',0);

d = ['74-108mm-delay--6ps\results-maxQuad-' num2str(maxQuad) '-Res-' ...
    num2str(Res) '-maxX-' num2str(maxX) '-Xstep-' num2str(XStep) '-phiStep-' num2str(phiStep)];
d = ['196-291,000mm-delay-1214ps\results-maxQuad-' num2str(maxQuad) '-Res-' ...
    num2str(Res) '-maxX-' num2str(maxX) '-Xstep-' num2str(XStep) '-phiStep-' num2str(phiStep)];


if ~exist(d,'dir')
    mkdir(d)
end

for Rvec = 0.1:0.1:2
    
%     directory = 'C:\Users\Public\Documents\archived-data\';
%     dispstat('makegrid','init','timestamp','keepthis',0);
%     mainCalcGridsAndPatterns('Directory',directory,'R',Rvec,'maxQuad',maxQuad,'Resolution',Res,'maxX',maxX,'XStep',XStep,'PhiStep',phiStep);
    directory = ['C:\Users\Public\Documents\archived-data\PFunction-R-' num2str(Rvec) ...
        '-maxQuad-' num2str(maxQuad) '-Resolution-' num2str(Res) '-maxX-' num2str(maxX) '-Xstep-' num2str(XStep) '-phiStep-' num2str(phiStep)];
    [xGrid,phiGrid] = getGridsFromFilenames(directory);
    tic;
    dispstat('calculate P function','init','timestamp','keepthis',0);
    [P,sigmaP,QuadVals] = PFunctionFromData(selX2,selTheta,directory,xGrid,phiGrid);
    toc;
    save([d '\P-' num2str(Rvec) '-maxQuad-' num2str(maxQuad) '-Resolution-' num2str(Res) '.mat'],'P','sigmaP','QuadVals','xGrid','phiGrid','Rvec');
    plotWigner(P,'PQ',QuadVals,'ZString',['P(q,p) with R = ' num2str(Rvec)],'Filename',[d '\P-R-' num2str(Rvec)]);
    close();
    plotWigner(sigmaP,'PQ',QuadVals,'ZString',['\sigma P(q,p) with R = ' num2str(Rvec)],'Filename',[d '\sigmaP-R-' num2str(Rvec)]);
    close();
    
end

Contents = dir([d '/*.mat']);
name = {Contents.name};
[Pmax,sigmaPmax,Rvec,width,meanPhase,meanAmp,meanPhaseBinned,meanAmpBinned,varPhase,varAmp,varPhaseBinned,varAmpBinned,PhotonNr,PhotonNrBinned,g1,sigNeg] = ...
    deal(zeros(length(Contents),1));
for i = 1:length(Contents)
    filename = cell2mat(name(i));
    RToken = regexpi(filename,'P-([-0123456789.]*)-','tokens');
    Rvec(i) = str2double(cell2mat(RToken{1}));
    load([d '/' filename],'P','sigmaP','QuadVals');
    % make sure P is normalized  
    sumP = sum(sum(P));
    P=P./sumP;
    sigmaP=sigmaP./sumP;
    [meanPhas,meanAm,meanPhaseBinne,meanAmpBinne,varPhas,varAm,varPhaseBinne,varAmpBinne,PhotonN,PhotonNrBinne,g,sigNe] = ...
        ReturnMomentsFromP( P, sigmaP, QuadVals, 30,[d '/R-' num2str(Rvec(i))] );
    meanPhase(i) = meanPhas; meanAmp(i) = meanAm; meanPhaseBinned(i) = meanPhaseBinne;
    meanAmpBinned(i) = meanAmpBinne; varPhase(i) = varPhas; varAmp(i) = varAm; 
    varPhaseBinned(i) = varPhaseBinne; varAmpBinned(i) = varAmpBinne; PhotonNr(i) = PhotonN; 
    PhotonNrBinned(i) = PhotonNrBinne; g1(i) = g; 
    if isempty(sigNe)
        sigNeg(i) = 0;
    else
        sigNeg(i) = sigNe;
    end
    Pmax(i) = max(max(P));
    sigmaPmax(i) = max(max(sigmaP));
    prq = sum(P,2);
    prq = prq * max(max(P))/max(prq);
    try
        width(i) = fwhm(QuadVals,prq);
    catch
        width(i) = 0;
    end

end

figure(1);
semilogy(Rvec,Pmax,'o-');
xlabel('R');
hold on;
semilogy(Rvec,sigmaPmax,'o-');
semilogy(Rvec,sigmaPmax./Pmax,'o-');
semilogy(Rvec,width,'o-');
l=legend('max(P)','max(\sigma P)','max(\sigma P) / max(P)','FWHM(P)','location','bestoutside');
graphicsSettings;
l.FontSize = 9;
savefig([d '/R-comparison.fig']);
print([d '/R-comparison.png'],'-dpng','-r300');

figure(2);
xlabel('R');
plot(Rvec,meanPhase,Rvec,meanPhaseBinned,Rvec,varPhase,Rvec,varPhaseBinned,Rvec,...
    g1,Rvec,sigNeg,'-');
l = legend('meanPhase','meanPhaseBinned','varPhase',...
    'varPhaseBinned','g1','sigNeg','location','bestoutside');
graphicsSettings;
l.FontSize = 9;
savefig([d '/R-comparison-phase-g1-sigNeg.fig']);
print([d '/R-comparison-phase-g1-sigNeg.png'],'-dpng','-r300');

figure(3);
xlabel('R');
plot(Rvec,meanAmp,Rvec,meanAmpBinned,Rvec,varAmp,Rvec,varAmpBinned,Rvec,...
    PhotonNr,Rvec,PhotonNrBinned,'-');
l = legend('meanAmp','meanAmpBinned','varAmp','varAmpBinned','PhotonNr','PhotonNrBinned','location','bestoutside');
graphicsSettings;
l.FontSize = 9;
savefig([d '/R-comparison-Amplitude-PhotonNr.fig']);
print([d '/R-comparison-Amplitude-PhotonNr.png'],'-dpng','-r300');

save([d '\R-comparison-maxQuad-' num2str(maxQuad) '-Res-' num2str(Res) '.mat'],...
    'Rvec','Pmax','sigmaPmax','width','Rvec','meanPhase','meanPhaseBinned','varPhase',...
    'varPhaseBinned','g1','sigNeg','meanAmp','meanAmpBinned','varAmp','varAmpBinned','PhotonNr','PhotonNrBinned');


