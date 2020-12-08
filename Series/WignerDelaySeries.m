function [] = WignerDelaySeries(varargin)
%get Wigner function for one delay, several postselection parameters, from
%data located in 'post-data'. 
%% Validate and parse input arguments
p = inputParser;
% directory where the Wigner tables are stored
defaultDirectory = 'C:\Users\Carolin Lüders\Documents\archived-data\Wigner';
%'C:\Users\lab\Documents\@archived-data\Wigner';
addParameter(p,'Directory',defaultDirectory,@isstr);
% maximum Fock state for density matrix
defaultMaxFockState = 60;
addParameter(p,'MaxFockState',defaultMaxFockState,@isnumeric);
% Iterations for density matrix
defaultIterations = 200;
addParameter(p,'Iterations',defaultIterations,@isnumeric);
% The position in mm of the zero delay
defaultZeroDelay = 108.96;
addParameter(p,'ZeroDelay',defaultZeroDelay,@isnumeric);
% Whether already existent rho and WF are used 
defaultLoadExistent = false;
addParameter(p,'LoadExistent',defaultLoadExistent,@islogical);
% Whether the WF is smoothed after calculating it  
defaultSmooth = false;
addParameter(p,'Smooth',defaultSmooth,@islogical);
% Half width for the moving average smoothing
defaultSmoothWidth = 2;
addParameter(p,'SmoothWidth',defaultSmoothWidth,@isnumeric);
defaultSelParams = struct('Type','fullcircle','Position',[7,0.4]);
addParameter(p,'SelectionParameters',defaultSelParams,@isstruct);
defaultRemoveModulation = false;
addParameter(p,'RemoveModulation',defaultRemoveModulation,@islogical);
defaultVaryAPS = false;
addParameter(p,'VaryAPS',defaultVaryAPS,@islogical);
defaultRange = [0 20];
addParameter(p,'Range',defaultRange,@isvector);
parse(p,varargin{:});
c = struct2cell(p.Results);
[directory,Iterations,loadExistent,maxFockState,range,remMod,selParams,smooth,smoothWidth,varyAPS,zeroDelay] = c{:};


%% preparation
selStr = selParamsToStr(selParams);
foldername = ['Wignerplots-',selStr,'-remMod-',...
            num2str(remMod),'-range-',num2str(range),'-varyAPS-',num2str(varyAPS)];
if ~exist([pwd foldername],'dir')
    mkdir(foldername)
end
filestruct = dir('mat-data\*.mat');
files = {filestruct.name};
filenamePlot = [foldername '\Iterations-' num2str(Iterations) '-maxN-' num2str(maxFockState) '-smooth-' num2str(smooth)];
MaxQuad = 20;
Resolution = 0.125;
[Q,P,varQ,varP,n,Delay,meanPhases,meanAmps,varPhases,varAmps,meanAbsPhases] = deal(zeros(length(files),1));

%% Iterate through data files
for i = 1:length(files)
    %% Load data
    C = strsplit(files{i},'.');
    filename = C{1};
    postFilename =  [filename,'-',selStr,'-remMod-',...
            num2str(remMod),'-range-',num2str(range),'-varyAPS-',num2str(varyAPS),'.mat'];
    filenameFig = [foldername '\' postFilename '-Iterations-' num2str(Iterations) '-maxN-' num2str(maxFockState) '-smooth-' num2str(smooth)];
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
    
    dispstat(['load file ' filename],'timestamp',0);
    if loadExistent
        load(['post-data/',postFilename],'WF','rho'); 
    else
        load(['post-data/',postFilename],'selX','selTheta');
        dispstat('compute rho','timestamp',0);
        rho = computeDensityMatrix(selX,selTheta,'MaxFockState',maxFockState,'Iterations',Iterations);
        dispstat('compute WF','timestamp',0);
        tic;
        WF = mainWignerFromRho(rho,'Directory',directory);
        toc;
        save(['post-data\' postFilename],'rho','WF','-append');
    end
       
    if smooth
        WF = moving_average(WF,smoothWidth,1);
        WF = moving_average(WF,smoothWidth,2);
    end
    plotWigner(WF,'Style','2D','Colormap','hot');
    title(['Delay =' num2str(delay) 'ps,' char(10) num2str(delayMm) 'mm']);       
    savefig([filenameFig '-WF.fig']);
    print([filenameFig '-WF.png'],'-dpng');
    close all;
    plotRho(rho);
    title(['Delay =' num2str(delay) 'ps, ' num2str(delayMm) 'mm']);
    savefig([filenameFig '-rho.fig']);
    print([filenameFig '-rho.png'],'-dpng');
    close all;
    [meanPhase,meanAmp,~,~,varPhase,varAmp,~,~,meanAbsPhase] = ReturnPhaseAndAmplitudeWigner( real(WF),...
    MaxQuad, Resolution,100,filenameFig );
    meanPhases(i) = meanPhase;
    meanAmps(i) = meanAmp;
    varPhases(i) = varPhase;
    varAmps(i) = varAmp;
    meanAbsPhases(i) = meanAbsPhase;
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
meanAbsPhases = meanAbsPhases(I);
save([foldername '\Wignerresults.mat'],'Delay','Q','P','varQ','varP','n','meanPhases','meanAmps','varPhases','varAmps');

plot(Delay,Q,'-o',Delay,P,'-o',Delay,meanAmps,'-o',Delay,meanPhases,'-o',Delay,meanAbsPhases,'-o');
legend('<Q>','<P>','<r>','<\phi>','<|\phi|>','location','bestoutside');
xlabel('Delay (ps)');
ylabel('Mean values');
graphicsSettings;
savefig([filenamePlot '-Amplitudes.fig']);
print([filenamePlot '-Amplitudes.png'],'-dpng');
ylim([-1 5]);
close all;

plot(Delay,varQ,'-o',Delay,varP,'-o',Delay,varAmps,'-o',Delay,varPhases,'-o');
legend('Var(Q)','Var(P)','Var(r)','Var(\phi)','location','bestoutside');
xlabel('Delay (ps)');
ylabel('Variance');
graphicsSettings;
savefig([filenamePlot '-Variance.fig']);
print([filenamePlot '-Variance.png'],'-dpng');
close all;

plot(Delay,varAmps,'-o');
legend('Var(r)','location','bestoutside');
xlabel('Delay (ps)');
ylabel('Variance');
graphicsSettings;
savefig([filenamePlot '-VarianceR.fig']);
print([filenamePlot '-VarianceR.png'],'-dpng');
close all;

plot(Delay,n,'-o');
legend('n','location','bestoutside');
xlabel('Delay (ps)');
ylabel('Photon Number');
ylim([5 12]);
graphicsSettings;
savefig([filenamePlot '-nPhotons.fig']);
print([filenamePlot '-nPhotons.png'],'-dpng');
close all;

end

