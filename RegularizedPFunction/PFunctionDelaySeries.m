function [] = PFunctionDelaySeries(varargin)
%compute regularized P functions for a delay series, several postselection parameters, from
%data located in 'post-data', and compute expectation values. 
%% Validate and parse input arguments
p = inputParser;
% The position in mm of the zero delay
defaultZeroDelay = 108.96;
addParameter(p,'ZeroDelay',defaultZeroDelay,@isnumeric);
% Whether already existent rho and WF are used 
defaultLoadExistent = false;
addParameter(p,'LoadExistent',defaultLoadExistent,@islogical);
defaultSelParams = struct('Type','fullcircle','Position',[7,0.4]);
addParameter(p,'SelectionParameters',defaultSelParams,@isstruct);
defaultRemoveModulation = false;
addParameter(p,'RemoveModulation',defaultRemoveModulation,@islogical);
defaultVaryAPS = false;
addParameter(p,'VaryAPS',defaultVaryAPS,@islogical);
defaultRange = [0 20];
addParameter(p,'Range',defaultRange,@isvector);
defaultPlot = true;
addParameter(p,'Plot',defaultPlot,@islogical);
%parameters used for the pattern function for the P function 
defaultMaxQuad = 20; %max abs value of the new quadrature coordinates 
addParameter(p,'MaxQuad',defaultMaxQuad,@isnumeric);
defaultMaxX= 20; %max abs value of the old quadrature coordinates 
addParameter(p,'MaxX',defaultMaxX,@isnumeric);
defaultPhiStep= 0.1; %step size of the phase grid
addParameter(p,'PhiStep',defaultPhiStep,@isnumeric);
defaultRvalue= 0.7; %filter Parameter R
addParameter(p,'Rvalue',defaultRvalue,@isnumeric);
defaultXStep= 1; %step size of the old quadrature coordinates grid
addParameter(p,'XStep',defaultXStep,@isnumeric);
defaultResolution= 1; %step size of the new quadrature coordinates grid
addParameter(p,'Resolution',defaultResolution,@isnumeric);
defaultNorm= 1; %normalization of the vacuum standard deviation of quadratures. Is 1 for P functions
addParameter(p,'Norm',defaultNorm,@isnumeric);
defaultPatternDir = 'C:\Users\Public\Documents\archived-data'; %directory with the pattern functions
% 'C:\Users\lab\Documents\@archived-data'
addParameter(p,'PatternDir',defaultPatternDir,@isstr);
parse(p,varargin{:});
c = struct2cell(p.Results);
[loadExistent,maxQuad,maxX,norm,patternDir,phiStep,plotOption,range,remMod,res,rvalue,selParams,varyAPS,XStep,zeroDelay] = c{:};

%directory with the pattern functions
directory = [patternDir '\PFunction-R-' num2str(rvalue) ...
        '-maxQuad-' num2str(maxQuad) '-Resolution-' num2str(res) '-maxX-' num2str(maxX) '-Xstep-' num2str(XStep) '-phiStep-' num2str(phiStep)];
[xGrid,phiGrid] = getGridsFromFilenames(directory);

%% preparation
selStr = selParamsToStr(selParams);
foldername = ['Pfunctionplots-',selStr,'-remMod-',...
            num2str(remMod),'-range-',num2str(range),'-varyAPS-',num2str(varyAPS),'-R-' num2str(rvalue) ...
        '-maxQuad-' num2str(maxQuad) '-Resolution-' num2str(res) '-maxX-' num2str(maxX) '-Xstep-' num2str(XStep) '-phiStep-' num2str(phiStep)];
if ~exist([pwd foldername],'dir')
    mkdir(foldername)
end

% Get the unique filename beginnings from post-data
postfilestruct = dir('post-data\*.mat');
postfiles = {postfilestruct.name};
postfilenames = cell(size(postfiles));
for i = 1:length(postfiles)
    C = strsplit(postfiles{i},'-type');
    postfilenames{i} = C{1};
end
files = unique(postfilenames);

[Delay,DelayMm,Pmax,sigmaPmax,meanPhase,meanAmp,meanPhaseBinned,meanAmpBinned,...
    varPhase,circVar1,circVar2,circVar1Err,circVar2Err,varAmp,varPhaseBinned,varAmpBinned,PhotonNr,PhotonNrBinned,g1,sigNeg,...
    phaseErr,ampErr,varPhaseErr,varAmpErr,PhotonNrErr,maxQ] = deal(zeros(length(files),1));

%% Iterate through data files
for i = 1:length(files)
    %% Load data
    postFilename =  [files{i},'-',selStr,'-remMod-',...
    num2str(remMod),'-range-',num2str(range),'-varyAPS-',num2str(varyAPS),'.mat'];
    filenameFig = [foldername '\' postFilename ];
    dispstat(['Loading ',files{i},' ...'],'timestamp',0);
    clear X1 X2 X3 theta piezoSign
    clear O1 O2 O3 oTheta selX selTheta;
    clear P sigmaP;
   
    % get delay from the file name with format xx-yymm, where xx is the
    % filenumber and yy is the delay which can also be negative and start
    % with a minus sign
    delayToken = regexpi(files{i},'([-0123456789,-]*)mm','tokens');
    delay = cell2mat(delayToken{1});
    numberToken = regexpi(delay,'^([0123456789,]*)-','tokens'); 
    number = cell2mat(numberToken{1});
    delay = strrep(delay,[number '-'],'');
    delay = strrep(delay,',','.');
    delay = str2double(delay);
    DelayMm(i) = delay;
    c = 299792458; % in m/s
    delay = 2*(delay-zeroDelay)/1000/c*10^12; %delay in ps
    Delay(i) = delay;
    
    dispstat(['load file ' files{i}],'timestamp',0);
    if loadExistent
        load(['post-data/',postFilename],'P','sigmaP','QuadVals'); 
    else
        load(['post-data/',postFilename],'selX','selTheta');
        %Normalization: If the standard deviation of the vacuum is set to 1 
        if norm == 1
            selX2 = sqrt(2)*selX;
        end
        dispstat('compute P function','timestamp',0);
        tic;
        dispstat('calculate P function','init','timestamp','keepthis',0);
        [P,sigmaP,QuadVals] = PFunctionFromData(selX2,selTheta,directory,xGrid,phiGrid);
        toc;
        save(['post-data\' postFilename],'P','sigmaP','QuadVals','-append');
    end
    
    % make sure P is normalized  
    sumP = sum(sum(P));
    P=P./sumP;
    sigmaP=sigmaP./sumP;
    
    if plotOption       
        plotWigner(P,'PQ',QuadVals,'ZString','P_{\Omega}(q,p)','Filename',[filenameFig '-P-']);
        close();
        plotWigner(sigmaP,'PQ',QuadVals,'ZString','\sigma P_{\Omega}(q,p)','Filename',[filenameFig '-sigmaP-']);
        close();
    end
    
    [meanPhas,meanAm,meanPhaseBinne,meanAmpBinne,varPhas,circVa1,circVa2,...
        circVar1Er,circVar2Er,varAm,varPhaseBinne,...
        varAmpBinne,PhotonN,PhotonNrBinne,g,sigNe,phaseEr,ampEr,varPhaseEr,varAmpEr,PhotonNrEr,maxq] = ...
        ReturnMomentsFromP( P, sigmaP, QuadVals, 30,filenameFig,'Plot',plotOption );
    meanPhase(i) = meanPhas; meanAmp(i) = meanAm; meanPhaseBinned(i) = meanPhaseBinne;
    meanAmpBinned(i) = meanAmpBinne; varPhase(i) = varPhas;
    circVar1(i) = circVa1; circVar2(i) = circVa2; ...
    circVar1Err(i) = circVar1Er; circVar2Err(i) = circVar2Er; varAmp(i) = varAm; 
    varPhaseBinned(i) = varPhaseBinne; varAmpBinned(i) = varAmpBinne; PhotonNr(i) = PhotonN; 
    PhotonNrBinned(i) = PhotonNrBinne; g1(i) = g; phaseErr(i) = phaseEr; ampErr(i) = ampEr; ...
        varPhaseErr(i) = varPhaseEr; varAmpErr(i) = varAmpEr; PhotonNrErr(i) = PhotonNrEr;
    maxQ(i) = maxq;
    if isempty(sigNe)
        sigNeg(i) = 0;
    else
        sigNeg(i) = sigNe;
    end
    Pmax(i) = max(max(P));
    sigmaPmax(i) = max(max(sigmaP));
end

[Delay,I]=sort(Delay);
DelayMm = DelayMm(I);
meanPhase = meanPhase(I); meanAmp = meanAmp(I); meanPhaseBinned = meanPhaseBinned(I);
meanAmpBinned = meanAmpBinned(I); varPhase = varPhase(I); 
circVar1 = circVar1(I); circVar2 = circVar2(I); 
circVar1Err = circVar1Err(I); circVar2Err = circVar2Err(I); varAmp = varAmp(I); 
varPhaseBinned = varPhaseBinned(I); varAmpBinned = varAmpBinned(I); PhotonNr = PhotonNr(I); 
PhotonNrBinned = PhotonNrBinned(I); g1 = g1(I);
Pmax = Pmax(I); sigmaPmax = sigmaPmax(I); sigNeg = sigNeg(I);
phaseErr = phaseErr(I); ampErr = ampErr(I); varPhaseErr = varPhaseErr(I); ...
varAmpErr = varAmpErr(I); PhotonNrErr = PhotonNrErr(I); maxQ = maxQ(I);
save([foldername '\Pfunctionresults.mat'],'Delay','DelayMm','Pmax','sigmaPmax','meanPhase','meanPhaseBinned','varPhase',...
    'varPhaseBinned','g1','sigNeg','meanAmp','meanAmpBinned','circVar1','circVar2','varAmp','varAmpBinned','PhotonNr','PhotonNrBinned',...
    'phaseErr','ampErr','varPhaseErr','varAmpErr','PhotonNrErr','circVar1Err','circVar2Err','maxQ');

figure(1);
plot(Delay,Pmax,'o-');
xlabel('Delay (ps)');
hold on;
plot(Delay,sigmaPmax,'o-');
plot(Delay,sigmaPmax./Pmax,'o-');
l=legend('max(P)','max(\sigma P)','max(\sigma P) / max(P)','location','bestoutside');
graphicsSettings;
l.FontSize = 9;
savefig([foldername '/Delay-Pmax-sigmaP.fig']);
print([foldername '/Delay-Pmax-sigmaP.png'],'-dpng','-r300');
close();

figure(2);
xlabel('R');
plot(Delay,meanPhase,Delay,meanPhaseBinned,Delay,varPhase,Delay,varPhaseBinned,Delay,...
    g1,Delay,sigNeg,Delay,circVar1,Delay,circVar2,'o-');
l = legend('meanPhase','meanPhaseBinned','varPhase',...
    'varPhaseBinned','g1','sigNeg','circVar1','circVar2','location','bestoutside');
graphicsSettings;
l.FontSize = 9;
savefig([foldername '/Delay-phase-g1-sigNeg.fig']);
print([foldername '/Delay-comparison-phase-g1-sigNeg.png'],'-dpng','-r300');
close();

figure(3);
xlabel('R');
plot(Delay,meanAmp,Delay,meanAmpBinned,Delay,varAmp,Delay,varAmpBinned,Delay,...
    PhotonNr,Delay,PhotonNrBinned,'o-');
l = legend('meanAmp','meanAmpBinned','varAmp','varAmpBinned','PhotonNr','PhotonNrBinned','location','bestoutside');
graphicsSettings;
l.FontSize = 9;
savefig([foldername '/Delay-Amplitude-PhotonNr.fig']);
print([foldername '/Delay-Amplitude-PhotonNr.png'],'-dpng','-r300');
close();



end

