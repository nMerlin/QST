%testing P function for different Rs
%Normalization: If the standard deviation of the vacuum is set to 1 
selX2 = sqrt(2)*selX;
maxQuad = 20;
Res = 1;

%d = '74-108mm-delay--6ps';
d = '196-291,000mm-delay-1214ps';

for Rvec = 0.1:0.1:2
    
%     directory = 'C:\Users\Public\Documents\archived-data\';
%     mainCalcGridsAndPatterns('Directory',directory,'R',Rvec);
    directory = ['C:\Users\Public\Documents\archived-data\PFunction-R-' num2str(Rvec) '-maxQuad-' num2str(maxQuad) '-Resolution-' num2str(Res)'];
    [xGrid,phiGrid] = getGridsFromFilenames(directory);
    [P,sigmaP,QuadVals] = PFunctionFromData(selX2,selTheta,directory,xGrid,phiGrid);
    save([d '\P-' num2str(Rvec) '-maxQuad-' num2str(maxQuad) '-Resolution-' num2str(Res) '.mat'],'P','sigmaP','QuadVals','xGrid','phiGrid','Rvec');
    plotWigner(P,'PQ',QuadVals,'ZString',['P(q,p) with R = ' num2str(Rvec)],'Filename',[d '\P-R-' num2str(Rvec)]);
    close();
    plotWigner(sigmaP,'PQ',QuadVals,'ZString',['\sigma P(q,p) with R = ' num2str(Rvec)],'Filename',[d '\sigmaP-R-' num2str(Rvec)]);
    close();
    
end

Contents = dir([d '/*.mat']);
name = {Contents.name};
[Pmax,sigmaPmax,Rvec,width] = deal(zeros(length(Contents),1));
for i = 1:length(Contents)
    filename = cell2mat(name(i));
    RToken = regexpi(filename,'P-([-0123456789.]*)-','tokens');
    Rvec(i) = str2double(cell2mat(RToken{1}));
    load([d '/' filename],'P','sigmaP','QuadVals');
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

semilogy(Rvec,Pmax,'o-');
xlabel('R');
hold on;
semilogy(Rvec,sigmaPmax,'o-');
semilogy(Rvec,sigmaPmax./Pmax,'o-');
semilogy(Rvec,width,'o-');
legend('max(P)','max(\sigma P)','max(\sigma P) / max(P)','FWHM(P)','location','best');
graphicsSettings;
savefig([d '/R-comparison.fig']);
print([d '/R-comparison.png'],'-dpng','-r300');
