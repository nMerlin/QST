%get Wigner function for one delay, several postselection parameters
directory = 'C:\Users\Carolin Lüders\Documents\archived-data\Wigner';
zeroDelay = 108.96;

if ~exist([pwd 'Wignerplots'],'dir')
    mkdir('Wignerplots')
end

foldername = 'Wignerplots';
filestruct = dir('post-data\*.mat');
files = {filestruct.name};

MaxQuad = 20;
Resolution = 0.125;
[Q,P,varQ,varP,n,Delay,meanPhases,meanAmps,varPhases,varAmps] = deal(zeros(length(files),1));

%% Iterate through data files

for i = 1:length(files)
    %% Load data
    C = strsplit(files{i},'.');
    filename = C{1};
    dispstat(['Loading ',files{i},' ...'],'timestamp',0);
    clear X1 X2 X3 theta piezoSign
    clear O1 O2 O3 oTheta selX selTheta;
    clear rho WF;
   
    % get delay from the file name with format xx-yymm, where xx is the
    % filenumber and yy is the delay which can also be negative and start
    % with a minus sign
    delayToken = regexpi(filename,'([-0123456789,-]*)mm','tokens');
    delay = cell2mat(delayToken{1});
    numberToken = regexpi(delay,'^([0123456789,]*)-','tokens'); 
    number = cell2mat(numberToken{1});
    delay = strrep(delay,[number '-'],'');
    delay = strrep(delay,',','.');
    delay = str2double(delay);
    delayMm = delay;
    c = 299792458; % in m/s
    delay = 2*(delay-zeroDelay)/1000/c*10^12; %delay in ps
    Delay(i) = delay;
     
%     load(['post-data/',files{i}],'selX','selTheta');
%     
%         tic;
%         rho = computeDensityMatrix(selX,selTheta,'MaxFockState',60,'Iterations',200);
%         WF = mainWignerFromRho(rho,'Directory',directory);
%         toc;
%        % save(['post-data\' filename selStr '-remMod-0-range-0.3-varyAPS-0.mat'],'rho','WF','-append');
%        save(['post-data\' files{i}],'rho','WF','-append');
         load(['post-data/',files{i}],'WF','rho');  
        plotWigner(WF,'Style','2D','Colormap','hot');
        title(['Delay =' num2str(delay) 'ps, ' num2str(delayMm) 'mm']);
        savefig([foldername '\' filename '-WF.fig']);
         print([foldername '\' filename '-WF.png'],'-dpng');
         close all;
         plotRho(rho);
        title(['Delay =' num2str(delay) 'ps, ' num2str(delayMm) 'mm']);
        savefig([foldername '\' filename '-rho.fig']);
         print([foldername '\' filename '-rho.png'],'-dpng');
         close all;
         [meanPhase,meanAmp,~,~,varPhase,varAmp,~,~ ] = ReturnPhaseAndAmplitudeWigner( real(WF),...
             MaxQuad, Resolution,100,[foldername '\' filename] );
         meanPhases(i) = meanPhase;
         meanAmps(i) = meanAmp;
         varPhases(i) = varPhase;
         varAmps(i) = varAmp;
         [Qx,Qy,FullQx2,FullQy2,VarQx,VarQy,PhotonNr ] = ReturnQuadsWigner( real(WF), MaxQuad, Resolution );
       Q(i) = Qx;
       P(i) = Qy;
       varQ(i) = VarQx;
        varP(i) = VarQy;
         n(i) = PhotonNr; 
end

[Delay,I]=sort(Delay);
Q=Q(I);
P=P(I);
varQ=varQ(I);
varP=varP(I);
n=n(I);
meanPhases = meanPhases(I);
meanAmps = meanAmps(I);
varPhases = varPhases(I);
varAmps = varAmps(I);
save([foldername '\Wignerresults.mat'],'Delay','Q','P','varQ','varP','n','meanPhases','meanAmps','varPhases','varAmps');

plot(Delay,Q,'-o',Delay,P,'-o',Delay,meanAmps,'-o');
legend('<Q>','<P>','<r>','location','best');
xlabel('Delay (ps)');
ylabel('Amplitude');
graphicsSettings;
savefig([foldername '\Amplitudes.fig']);
print([foldername '\Amplitudes.png'],'-dpng');
close all;

plot(Delay,varQ,'-o',Delay,varP,'-o',Delay,varAmps,'-o',Delay,varPhases,'-o');
legend('Var_{Q}','Var_{P}','Var_{r}','Var(\phi)','location','best');
xlabel('Delay (ps)');
ylabel('Variance');
graphicsSettings;
savefig([foldername '\Variance.fig']);
print([foldername '\Variance.png'],'-dpng');
close all;

plot(Delay,n,'-o');
legend('n','location','best');
xlabel('Delay (ps)');
ylabel('Photon Number');
graphicsSettings;
savefig([foldername '\nPhotons.fig']);
print([foldername '\nPhotons.png'],'-dpng');
close all;

plot(Delay,meanPhases,'-o');
legend('<\phi>','location','best');
xlabel('Delay (ps)');
ylabel('Mean Phase');
graphicsSettings;
savefig([foldername '\meanPhase.fig']);
print([foldername '\meanPhase.png'],'-dpng');
close all;


