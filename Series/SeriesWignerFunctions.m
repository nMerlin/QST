%get Wigner function for one delay, several postselection parameters
t = 0.4;
listOfParams = struct('Type',{'fullcircle'},'Position',{[0.5 t],[1,t],[1.5 t],...
[2 t],[3 t],[4,t],[5 t],[6 t],[7 t],[8 t],[9 t],[10 t],[11 t],[12 t],[13 t]});
% filename = '278--4,600mm-493,904mW-4mW-LOwithDL-';
filename = '278--4,600mm-493,904mW-4mW-LOwithDL-corrRemove-orthwidth-0.02-Wigner';
%foldername = ['Wignerfunctions-278--4.6mm-radiusSeries-thickness-' num2str(t)];
foldername = ['Wignerfunctions-278--4.6mm-radiusSeries-thickness-' num2str(t) '-corrRemove-orthwidth-0.02'];
if ~exist([pwd foldername],'dir')
    mkdir(foldername)
end

for iParams = 1:length(listOfParams)
        selParams = listOfParams(iParams);
         selStr = selParamsToStr(selParams);
%         load(['post-data\' filename selStr '-remMod-0-range-0.3-varyAPS-0.mat']);
        [selX,selTheta] = selectRegion(O1,O2,O3,oTheta,selParams);
        tic;
        rho = computeDensityMatrix(selX,selTheta,'MaxFockState',60);
        WF = mainWignerFromRho(rho);
        toc;
       % save(['post-data\' filename selStr '-remMod-0-range-0.3-varyAPS-0.mat'],'rho','WF','-append');
       save([foldername '\' filename selStr '-remMod-0-range-0.3-varyAPS-0.mat'],'rho','WF');
        plotWigner(WF,'Style','2D');
        savefig([foldername '\' filename selStr '-remMod-0-range-0.3-varyAPS-0-WF.fig']);
         print([foldername '\' filename selStr '-remMod-0-range-0.3-varyAPS-0-WF.png'],'-dpng');
         close all;
end

MaxQuad = 20;
Resolution = 0.125;
[Q,P,varQ,varP,n,Aps] = deal(zeros(length(listOfParams)));

for iParams = 1:length(listOfParams)
        selParams = listOfParams(iParams);
         selStr = selParamsToStr(selParams);
         Aps(iParams) = selParams.Position(1);
        %load(['post-data\' filename selStr '-remMod-0-range-0.3-varyAPS-0.mat'],'WF'); 
        load([foldername '\' filename selStr '-remMod-0-range-0.3-varyAPS-0.mat'],'WF'); 
       [Qx,Qy,FullQx2,FullQy2,VarQx,VarQy,PhotonNr ] = ReturnQuadsWigner( WF, MaxQuad, Resolution );
       Q(iParams) = Qx;
       P(iParams) = Qy;
       varQ(iParams) = VarQx;
        varP(iParams) = VarQy;
         n(iParams) = PhotonNr;     
end

plot(Aps,Q,'-o',Aps,P,'-o');
legend('<Q>','<P>','location','best');
xlabel('A_{PS}');
ylabel('Amplitude');
graphicsSettings;
savefig([foldername '\Amplitudes.fig']);
print([foldername '\Amplitudes.png'],'-dpng');
close all;

plot(Aps,varQ,'-o',Aps,varP,'-o');
legend('Var_{Q}','Var_{P}','location','best');
xlabel('A_{PS}');
ylabel('Variance');
graphicsSettings;
savefig([foldername '\Variance.fig']);
print([foldername '\Variance.png'],'-dpng');
close all;

plot(Aps,n,'-o');
legend('n','location','best');
xlabel('A_{PS}');
ylabel('Photon Number');
graphicsSettings;
savefig([foldername '\nPhotons.fig']);
print([foldername '\nPhotons.png'],'-dpng');
close all;


