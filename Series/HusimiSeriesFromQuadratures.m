function HusimiSeriesFromQuadratures( chAssign, varargin)
%computes and plots Husimi functions from quadrature data located in the folder 'mat-data'. 
%'ChannelAssignment': This sets which channel is target, 
%       which is the postselection channel that is piezo modulated fast for 
%       the phase computation and which is the last postselection channel, 
%       according to [target,ps_piezo_fast,ps_piezo_slow]
% For Husimi funciton, only ps_piezo_fast and ps_piezo_slow are used. 
% Optional input arguments
%% Validate and parse input arguments
p = inputParser;
defaultParameter = 'Power'; % which parameter was changed during the series
addParameter(p,'Parameter',defaultParameter);
defaultSaveHusimi = false; % if the H, binsO1... are saved
addParameter(p,'SaveHusimi',defaultSaveHusimi,@islogical);
defaultLoadExistent = false; % if data is already there and can be loaded.
addParameter(p,'LoadExistent',defaultLoadExistent,@islogical);
defaultXUnit = 'mW'; % unit  of the parameter
addParameter(p,'XUnit',defaultXUnit);
defaultRange = 0.5; % if maxN-minN > range*medN, there are two fits made 
%for higher and for lower photon numbers.  
addParameter(p,'range',defaultRange);
defaultPlotErrorbars = false; % if true, errorbars are plotted 
addParameter(p,'PlotErrorbars',defaultPlotErrorbars,@islogical);
defaultLOpower = 4; % important to read out filename
addParameter(p,'LOpower',defaultLOpower,@isnumeric);
defaultShowLegend = true;
addParameter(p,'ShowLegend',defaultShowLegend,@islogical);
defaultFitMethod = 'NLSQ-LAR'; %this is the most stable method 
addParameter(p,'FitMethod',defaultFitMethod,@isstr);
defaultScale = true; %scales the O1 and O2 so they have the same photon number
addParameter(p,'Scale',defaultScale,@islogical);
defaultPlotOption = true;
addParameter(p,'PlotOption',defaultPlotOption,@islogical);
defaultMonteCarloError = true; %get errors from MonteCarlo error estimation
addParameter(p,'MonteCarloError',defaultMonteCarloError,@islogical);
defaultNMonteCarlo = 1000; % size of sample for MonteCarlo error estimation
addParameter(p,'NMonteCarlo',defaultNMonteCarlo,@isnumeric);
defaultPthr = 0; %Pthreshold (mW); if this is >0, the Power is plotted in units of PThr.
addParameter(p,'Pthr',defaultPthr,@isnumeric);
parse(p,varargin{:});
c = struct2cell(p.Results);
[fitMethod,loadExistent,LOpower,monteCarloError,NMonteCarlo,parameter,...
    plotErrorbars,plotOption,Pthr,range,saveHusimi,scale,showLegend,xUnit] = c{:};

%% Variables
dataStruct = struct; %will contain the quantities of interest 

%% Create data overview
dispstat('','init');
dispstat('Checking filenames ...','timestamp','keepthis');
Contents = dir('mat-data');
name = {Contents.name};

for iStruct =  1:length(Contents) 
    %get filename
    filename = cell2mat(name(iStruct));
    if strcmp(filename,'.') || strcmp(filename,'..') || strcmp(filename,'.txt')
        continue
    end
          
    dataStruct(iStruct).filename = filename;
    
    switch parameter
        case 'Power'    
            %currentToken = regexpi(filename,'([0123456789,]*)mW','tokens');
            currentToken = regexpi(filename,['([0123456789,]*)mW-' num2str(LOpower) 'mW'],'tokens');
             currentToken{1}=strrep(currentToken{1},',','.');
             dataStruct(iStruct).I = str2double(cell2mat(currentToken{1}));
        case 'delay'
            delayToken = regexpi(filename,'([-0123456789,-]*)mm','tokens');
            delay = cell2mat(delayToken{1});
            numberToken = regexpi(delay,'^([0123456789,]*)-','tokens'); 
            number = cell2mat(numberToken{1});
            delay = strrep(delay,[number '-'],'');
            delay = strrep(delay,',','.');
            delay = str2double(delay);
            c = 299792458; % in m/s
            delay = 2*delay/1000/c*10^12; %delay in ps   
            dataStruct(iStruct).I = delay;
        case 'no' 
            dataStruct(iStruct).I = 0;
        case 'Position'
            delayToken = regexpi(filename,'([-0123456789,-]*)mm','tokens');
            delay = cell2mat(delayToken{1});
            numberToken = regexpi(delay,'^([0123456789,]*)-','tokens'); 
             dataStruct(iStruct).number = str2double(cell2mat(numberToken{1}));
            number = cell2mat(numberToken{1});
            delay = strrep(delay,[number '-'],'');
            delay = strrep(delay,',','.');
            delay = str2double(delay);
            dataStruct(iStruct).I = delay;
    end
         
    if ~exist('Husimiplots','dir')
    mkdir('Husimiplots')
    end
    
    %%% Load data
    
    if loadExistent
         dispstat(['load already existent data from ' filename],...
        'timestamp','keepthis','notquiet'); 
        load(['mat-data\' filename],'O1','O2','iOrth','nPsFast','nPsSlow','nPsFastVec');
    else
       dispstat(['load all data from ' filename],...
        'timestamp','keepthis','notquiet'); 
        load(['mat-data\' filename]);
        if ~exist('nPsFast','var') 
            
             if ~isequal(size(X1,1),size(X2,1),size(X3,1))
                 X1 = X1(1:min([size(X1,1),size(X2,1),size(X3,1)]),:,:);
                 X2 = X2(1:min([size(X1,1),size(X2,1),size(X3,1)]),:,:);
                 X3 = X3(1:min([size(X1,1),size(X2,1),size(X3,1)]),:,:);
             end
            
            quadratures = zeros([size(X1) 3]);
            quadratures(:,:,:,1) = X1;
            quadratures(:,:,:,2) = X2;
            quadratures(:,:,:,3) = X3;
            XpsFast = quadratures(:,:,:,chAssign(2));
            XpsSlow = quadratures(:,:,:,chAssign(3));
            clear('quadratures');
            [~,nPsFast,nPsSlow] = nPhotons(XpsFast,XpsFast,XpsSlow);
            nPsFastVec = photonNumberVector(XpsFast);
            save(['mat-data\' filename],'nPsFast','nPsSlow','nPsFastVec','-append');
        end
        
        if ~exist('O1','var')  % make orthogonal quadratures 
            [a,b,c] = size(XpsFast);
            theta = reshape(XpsFast,[a*b c]); %not real theta, but needed for selectOrthogonal 
            [O1,O2,~,~,iOrth] = selectOrthogonal(XpsFast,XpsSlow,XpsSlow,theta,piezoSign,'width',0.05);
            save(['mat-data\' filename],'O1','O2','iOrth','-append');
        end
    end
    
    if scale && mean([nPsFast nPsSlow]) >=1
        % scaling works only properly for photon numbers >= 1
        O1 = O1*sqrt(mean([nPsFast nPsSlow])/nPsFast);
        O2 = O2*sqrt(mean([nPsFast nPsSlow])/nPsSlow);
        scaled = true;
    else 
        scaled = false; 
    end
  
    maxN = max(nPsFastVec(:));
    minN = min(nPsFastVec(:));
    medN = mean([maxN minN]);
    if maxN-minN < range*medN  || nPsFast <1     
    % If the fluctuation is small, there is probably only one state all the
    % time, so one fit is enough. 
        [H, binsO1, binsO2,r,nTherm,nThermErr,nCoherent,nCohErr,meanN,Coherence,CoherenceErr,poissonErrors,nRatio,g2] = ...
            plotHusimiAndCut(O1,O2,0,0,0,'Filename',['Husimiplots\' filename ...
            '-scaled-' num2str(scaled) '-showLegend-' num2str(showLegend)],...
            'ShowLegend',showLegend,'FitMethod',fitMethod,'Plot',plotOption,'MonteCarloError',false);
         close all; 
        if saveHusimi
            save(['mat-data\' filename],'H','binsO1','binsO2','-append');
        end
        rscMax = r;
        nRatioErr=0;
        g2Err = 0;
        meanNErr = 0;
        
        if monteCarloError
        % make Monte Carlo Error Propagation
            [nThermRand,nCoherentRand,CoherenceRand,nRatioRand,g2Rand,meanNRand] = deal(zeros(NMonteCarlo,1));
            fitFunction = fittype('0.5*WRes^2*(pi*(a1+1))^-1 *exp(-(x.^2 + b1)/(a1+1)) .* besseli(0,2*x*sqrt(b1)/(a1+1))','problem','WRes'); 
            parfor i = 1:NMonteCarlo
                Hrandom = normrnd(H,poissonErrors);               
                [~, ~, ~,~,nTh,~,nC,~,n,C,~,~,nR,g] = ...
                plotHusimiAndCut(O1,O2,Hrandom,binsO1,binsO2,'Filename','-',...
                'FitMethod',fitMethod,'Plot',false,'fitFunction',fitFunction,'MonteCarloError',true);
                nThermRand(i) = nTh;
                nCoherentRand(i) = nC;
                CoherenceRand(i) = C;
                nRatioRand(i) = nR;
                g2Rand(i) = g;
                meanNRand(i) = n;
            end
            nThermMean = mean(nThermRand);
            nThermErr = std(nThermRand);
            nCohMean = mean(nCoherentRand);
            nCohErr = std(nCoherentRand);
            CoherenceMean = mean(CoherenceRand);        
            CoherenceErr = std(CoherenceRand);
            nRatioMean = mean(nRatioRand);
            nRatioErr = std(nRatioRand);
            g2Err = std(g2Rand);
            meanNErr = std(meanNRand);
        end
        
        [nThermMax,nThermHigh,nThermLow] = deal(nTherm);
        [nCoherentMax,nCoherentHigh,nCoherentLow] = deal(nCoherent);
        [weightLow,weightHigh] = deal(0);
        [CoherenceMax,CoherenceHigh,CoherenceLow] = deal(Coherence);
        [nThermErrMax,nThermErrHigh,nThermErrLow] = deal(nThermErr);
        [nCohErrMax,nCohErrHigh,nCohErrLow] = deal(nCohErr);
        [CoherenceErrMax,CoherenceErrHigh,CoherenceErrLow] = deal(CoherenceErr);
        [nRatioLow,nRatioHigh]=deal(nRatio);
        [nRatioErrLow,nRatioErrHigh]=deal(nRatioErr);
        [g2Low,g2High]=deal(g2);
        [g2ErrLow,g2ErrHigh]=deal(g2Err);
        [meanNMax,meanNHigh,meanNLow] = deal(meanN);
        [meanNErrLow,meanNErrHigh]=deal(meanNErr);
    else
    % If the std is big, there is probably switching between states
    % over time, so we make a fit for the lower and the higher state
    % each. 
        nPsFastVec = nPsFastVec(iOrth);
        rangeLow = [0 medN];
        %rangeLow = [0 (medN+minN)/2];
        rangeHigh = [medN maxN];
        %rangeHigh = [(medN+maxN)/2 maxN];
        % for the low one 
        iSelLow = find(nPsFastVec >= min(rangeLow) & nPsFastVec <= max(rangeLow));
        O1s = O1(iSelLow);
        O2s = O2(iSelLow);
        [H, binsO1, binsO2,rLow,nThermLow,nThermErrLow,nCoherentLow,nCohErrLow,meanNLow,CoherenceLow,CoherenceErrLow,poissonErrors,nRatioLow,g2Low] = ...
            plotHusimiAndCut(O1s,O2s,0,0,0,'Filename',...
            ['Husimiplots\' filename '-scaled-' num2str(scaled) '-showLegend-' num2str(showLegend) '-LowPhotonNumber'],...
            'ShowLegend',showLegend,'FitMethod',fitMethod,'Plot',plotOption,'MonteCarloError',false);
         close all;
         nRatioErrLow=0;
         g2ErrLow = 0;
         meanNErrLow = 0;
         if monteCarloError
        % make Monte Carlo Error Propagation
            [nThermRand,nCoherentRand,CoherenceRand,nRatioRand,g2Rand,meanNRand] = deal(zeros(NMonteCarlo,1));
            fitFunction = fittype('0.5*WRes^2*(pi*(a1+1))^-1 *exp(-(x.^2 + b1)/(a1+1)) .* besseli(0,2*x*sqrt(b1)/(a1+1))','problem','WRes'); 
            parfor i = 1:NMonteCarlo
                Hrandom = normrnd(H,poissonErrors);               
                [~, ~, ~,~,nTh,~,nC,~,n,C,~,~,nR,g] = ...
                plotHusimiAndCut(O1,O2,Hrandom,binsO1,binsO2,'Filename','-',...
                'FitMethod',fitMethod,'Plot',false,'fitFunction',fitFunction,'MonteCarloError',true);
                nThermRand(i) = nTh;
                nCoherentRand(i) = nC;
                CoherenceRand(i) = C;
                nRatioRand(i) = nR;
                g2Rand(i) = g;
                meanNRand(i) = n;
            end
            nThermErrLow = std(nThermRand);
            nCohErrLow = std(nCoherentRand);    
            CoherenceErrLow = std(CoherenceRand);
            nRatioErrLow = std(nRatioRand);
            g2ErrLow = std(g2Rand);
            meanNErrLow = std(meanNRand);
        end
        
         % for the high one 
        iSelHigh = find(nPsFastVec >= min(rangeHigh) & nPsFastVec <= max(rangeHigh));
        O1s = O1(iSelHigh);
        O2s = O2(iSelHigh);
        [H, binsO1, binsO2, rHigh,nThermHigh,nThermErrHigh,nCoherentHigh,nCohErrHigh,meanNHigh,CoherenceHigh,CoherenceErrHigh,poissonErrors,nRatioHigh,g2High] = ...
            plotHusimiAndCut(O1s,O2s,0,0,0,'Filename',...
            ['Husimiplots\' filename '-scaled-' num2str(scaled) '-showLegend-' num2str(showLegend) '-HighPhotonNumber'],...
            'ShowLegend',showLegend,'FitMethod',fitMethod,'Plot',plotOption,'MonteCarloError',false);
         close all;
         nRatioErrHigh=0;
         g2ErrHigh = 0;
         meanNErrHigh = 0;
         if monteCarloError
        % make Monte Carlo Error Propagation
            [nThermRand,nCoherentRand,CoherenceRand,nRatioRand,g2Rand,meanNRand] = deal(zeros(NMonteCarlo,1));
            fitFunction = fittype('0.5*WRes^2*(pi*(a1+1))^-1 *exp(-(x.^2 + b1)/(a1+1)) .* besseli(0,2*x*sqrt(b1)/(a1+1))','problem','WRes'); 
            parfor i = 1:NMonteCarlo
                Hrandom = normrnd(H,poissonErrors);               
                [~, ~, ~,~,nTh,~,nC,~,n,C,~,~,nR,g] = ...
                plotHusimiAndCut(O1,O2,Hrandom,binsO1,binsO2,'Filename','-',...
                'FitMethod',fitMethod,'Plot',false,'fitFunction',fitFunction,'MonteCarloError',true);
                nThermRand(i) = nTh;
                nCoherentRand(i) = nC;
                CoherenceRand(i) = C;
                nRatioRand(i) = nR;
                g2Rand(i) = g;
                meanNRand(i) = n;
            end
            nThermErrHigh = std(nThermRand);
            nCohErrHigh = std(nCoherentRand);    
            CoherenceErrHigh = std(CoherenceRand);
            nRatioErrHigh = std(nRatioRand);
            g2ErrHigh = std(g2Rand);
            meanNErrHigh = std(meanNRand);
        end
        % their respective weights
        weightLow = length(iSelLow)/(length(iSelLow)+length(iSelHigh));
        weightHigh = length(iSelHigh)/(length(iSelLow)+length(iSelHigh));
        % weighted results
        r = weightLow*rLow + weightHigh*rHigh;
        nTherm = weightLow*nThermLow + weightHigh*nThermHigh;
        nCoherent = weightLow*nCoherentLow + weightHigh*nCoherentHigh;
        meanN = weightLow*meanNLow + weightHigh*meanNHigh;
        Coherence = weightLow*CoherenceLow + weightHigh*CoherenceHigh;
        nThermErr = weightLow*nThermErrLow + weightHigh*nThermErrHigh;
        nCohErr = weightLow*nCohErrLow + weightHigh*nCohErrHigh;
        CoherenceErr = weightLow*CoherenceErrLow + weightHigh*CoherenceErrHigh;
        % maximum results
        rscMax = max([rLow rHigh]);
        nThermMax = max([nThermLow nThermHigh]);
        nCoherentMax = max([nCoherentLow nCoherentHigh]);
        meanNMax = max([meanNLow meanNHigh]);
        CoherenceMax = max([CoherenceLow CoherenceHigh]);
        nThermErrMax = max([nThermErrHigh nThermErrLow]);
        nCohErrMax = max([nCohErrHigh nCohErrLow]);
        CoherenceErrMax = max([CoherenceErrHigh CoherenceErrLow]);
    end            
%     end  
    clear X1 X2 X3 theta piezoSign O1 O2 'H' 'binsO1' 'binsO2' 
    clear XpsFast XpsSlow nPsFast nPsSlow nPsFastVec
    dataStruct(iStruct).r = r;
    dataStruct(iStruct).nTherm = nTherm;
    dataStruct(iStruct).nCoherent = nCoherent;
    dataStruct(iStruct).meanN = meanN;  
    dataStruct(iStruct).Coherence = Coherence;
    dataStruct(iStruct).rMax = rscMax;
    dataStruct(iStruct).nThermMax = nThermMax;
    dataStruct(iStruct).nCoherentMax = nCoherentMax;
    dataStruct(iStruct).meanNMax = meanNMax; 
    dataStruct(iStruct).CoherenceMax = CoherenceMax;  
    dataStruct(iStruct).weightLow = weightLow;
    dataStruct(iStruct).weightHigh = weightHigh;
    dataStruct(iStruct).nThermHigh = nThermHigh;
    dataStruct(iStruct).nCoherentHigh = nCoherentHigh;
    dataStruct(iStruct).meanNHigh = meanNHigh; 
    dataStruct(iStruct).CoherenceHigh = CoherenceHigh;  
    dataStruct(iStruct).nThermLow = nThermLow;
    dataStruct(iStruct).nCoherentLow = nCoherentLow;
    dataStruct(iStruct).meanNLow = meanNLow;
    dataStruct(iStruct).CoherenceLow = CoherenceLow;     
    dataStruct(iStruct).nThermErr = nThermErr; 
    dataStruct(iStruct).nCohErr = nCohErr;  
    dataStruct(iStruct).CoherenceErr = CoherenceErr;
    dataStruct(iStruct).nThermErrMax = nThermErrMax;
    dataStruct(iStruct).nCohErrMax = nCohErrMax;
    dataStruct(iStruct).CoherenceErrMax = CoherenceErrMax;  
    dataStruct(iStruct).nRatioLow = nRatioLow; 
    dataStruct(iStruct).nRatioErrLow = nRatioErrLow; 
    dataStruct(iStruct).nRatioHigh = nRatioHigh; 
    dataStruct(iStruct).nRatioErrHigh = nRatioErrHigh; 
    dataStruct(iStruct).g2Low = g2Low; 
    dataStruct(iStruct).g2ErrLow = g2ErrLow; 
    dataStruct(iStruct).g2High = g2High; 
    dataStruct(iStruct).g2ErrHigh = g2ErrHigh;     
    dataStruct(iStruct).meanNErrLow = meanNErrLow;
    dataStruct(iStruct).meanNErrHigh = meanNErrHigh; 
end % iStruct

Is = cell2mat({dataStruct.I});
r = cell2mat({dataStruct.r});
nTherm = cell2mat({dataStruct.nTherm});
nCoherent = cell2mat({dataStruct.nCoherent});
meanN = cell2mat({dataStruct.meanN});
rMax = cell2mat({dataStruct.rMax});
nThermMax = cell2mat({dataStruct.nThermMax});
nCoherentMax = cell2mat({dataStruct.nCoherentMax});
meanNMax = cell2mat({dataStruct.meanNMax});
weightLow = cell2mat({dataStruct.weightLow});
weightHigh = cell2mat({dataStruct.weightHigh});
nThermHigh = cell2mat({dataStruct.nThermHigh});
nCoherentHigh = cell2mat({dataStruct.nCoherentHigh});
meanNHigh = cell2mat({dataStruct.meanNHigh});
nThermLow = cell2mat({dataStruct.nThermLow});
nCoherentLow = cell2mat({dataStruct.nCoherentLow});
meanNLow = cell2mat({dataStruct.meanNLow});
CoherenceLow = cell2mat({dataStruct.CoherenceLow});
CoherenceHigh = cell2mat({dataStruct.CoherenceHigh});
CoherenceMax = cell2mat({dataStruct.CoherenceMax});
Coherence = cell2mat({dataStruct.Coherence});
nThermErr = cell2mat({dataStruct.nThermErr});
nCohErr = cell2mat({dataStruct.nCohErr});
CoherenceErr = cell2mat({dataStruct.CoherenceErr});
nThermErrMax = cell2mat({dataStruct.nThermErrMax});
nCohErrMax = cell2mat({dataStruct.nCohErrMax});
CoherenceErrMax = cell2mat({dataStruct.CoherenceErrMax});
nRatioLow = cell2mat({dataStruct.nRatioLow});
nRatioErrLow = cell2mat({dataStruct.nRatioErrLow});
nRatioHigh = cell2mat({dataStruct.nRatioHigh});
nRatioErrHigh = cell2mat({dataStruct.nRatioErrHigh});
g2Low = cell2mat({dataStruct.g2Low});
g2ErrLow = cell2mat({dataStruct.g2ErrLow});
g2High = cell2mat({dataStruct.g2High});
g2ErrHigh = cell2mat({dataStruct.g2ErrHigh});
meanNErrLow = cell2mat({dataStruct.meanNErrLow});
meanNErrHigh = cell2mat({dataStruct.meanNErrHigh});

save(['Husimiplots\HusimiResultsScaled-FitMethod-' fitMethod '-scaled-' num2str(scaled) '-MCErr-' num2str(monteCarloError) '.mat'],'Is','r','nTherm','meanN','nCoherent',...
    'rMax','nThermMax','meanNMax','nCoherentMax','weightLow', 'weightHigh', 'nThermHigh', 'nCoherentHigh', 'meanNHigh', 'nThermLow', 'meanNLow','nCoherentLow',...
    'CoherenceHigh','CoherenceMax','Coherence','CoherenceLow','nThermErr',...
    'nCohErr','CoherenceErr','nThermErrMax','nCohErrMax','CoherenceErrMax',...
    'nRatioLow','nRatioErrLow','nRatioHigh','nRatioErrHigh',...
    'g2Low','g2ErrLow','g2High','g2ErrHigh','meanNErrLow','meanNErrHigh');
xlswrite(['Husimiplots\HusimiResultsScaled-' fitMethod '-scaled-' num2str(scaled) '-MCErr-' num2str(monteCarloError) '.xls'],[Is' r' nTherm' meanN' nCoherent' ...
    rMax' nThermMax' meanNMax' nCoherentMax' weightLow' weightHigh' nThermHigh' nCoherentHigh' meanNHigh' nThermLow' meanNLow' nCoherentLow' ...
     CoherenceHigh' CoherenceMax' Coherence' CoherenceLow' nThermErr' nCohErr',...
     CoherenceErr' nThermErrMax' nCohErrMax' CoherenceErrMax' ...
     nRatioLow', nRatioErrLow', nRatioHigh', nRatioErrHigh' ...
     g2Low', g2ErrLow', g2High', g2ErrHigh' meanNErrLow' meanNErrHigh']);

 %% plotting
if strcmp(parameter,'Power')
    parameter = 'Excitation power P';
end
fontsize = 17;
fontname = 'Arial';
msize = 8;
xmin = 3.5;
filenameOptions = ['-scaled-' num2str(scale) '-FitMethod-' fitMethod ...
    '-errorbars-' num2str(plotErrorbars) '-MCErr-' num2str(monteCarloError) '-Pthr-' num2str(Pthr)]; 
if Pthr > 0
    Is = Is/Pthr;
    xUnit = 'P_{thr}';
end

%% coherence
if plotErrorbars
    errorbar(Is,CoherenceLow,zeros(length(Is),1),CoherenceErr,'ok','linewidth',1.5,'markerSize',msize);
    f = gca;f.XScale = 'log'; hold on; 
    errorbar(Is,CoherenceHigh,zeros(length(Is),1),CoherenceErr,'*k','linewidth',1.5,'markerSize',msize);
else
    semilogx(Is,CoherenceLow,'ok',Is,CoherenceHigh,'*k','markerSize',10);
end
l=legend('$\mathcal{C}_\mathsf{low}$','$\mathcal{C}_\mathsf{high}$','location',...
    'northwest');
l.Interpreter = 'latex';
l.FontSize = fontsize+3;
ylabel('Quantum coherence');
xlabel([parameter ' (' xUnit ')']);
if Pthr > 0
    xlim([0.1 20]);
else 
    xlim([xmin max(Is)+10]);
end
graphicsSettings;
ax = gca;
set(ax,'FontSize',fontsize,'FontName',fontname);
savefig(['Husimiplots\CoherenceFromHusimiFunctions-semilogx' filenameOptions '.fig']);
print(['Husimiplots\CoherenceFromHusimiFunctions-semilogx' filenameOptions '.png'],'-dpng','-r300');
clf();

%% photon numbers all together
if plotErrorbars
    errorbar(Is,nThermLow,zeros(length(Is),1),nThermErr,'ok','linewidth',1.5,'markerSize',6);
    f = gca; f.XScale = 'log';f.YScale = 'log';hold on; 
    errorbar(Is,nCoherentLow,zeros(length(Is),1),nCohErr,'or','linewidth',1.5,'markerSize',6);
    errorbar(Is,nThermHigh,zeros(length(Is),1),nThermErr,'*k','linewidth',1.5,'markerSize',6);
    errorbar(Is,nCoherentHigh,zeros(length(Is),1),nCohErr,'*r','linewidth',1.5,'markerSize',6);    
    errorbar(Is,meanNLow,zeros(length(Is),1),meanNErrLow,':k','linewidth',1.5,'markerSize',10);
    errorbar(Is,meanNHigh,zeros(length(Is),1),meanNErrHigh,'-k','linewidth',1.5,'markerSize',10);
else
    loglog(Is,nThermLow,'ok',Is,nCoherentLow,'or',Is,nThermHigh,'*k',Is,nCoherentHigh,'*r','markerSize',8);
    hold on;
     loglog(Is,meanNLow,':k',Is,meanNHigh,'-k','linewidth',1.5);
end
% l = legend('$\bar n_\mathrm{low}$','$|\alpha_{0,\mathrm{\,low}}|^2$','$\bar n_\mathrm{high}$',...
%     '$|\alpha_{0,\mathrm{\,high}}|^2$','$n_{\mathrm{total,\,low}}$',...
%     '$n_{\mathrm{total,\,high}}$','location','northwest');
% l.FontSize = fontsize+3;
% sans-serif font: 
l = legend('$\bar \mathsf{n}_{\mathsf{low}}$','$\mathsf{|\alpha_{0,\,low}|^2}$','$\bar \mathsf{n}_{\mathsf{high}}$',...
    '$\mathsf{|\alpha_{0,\,high}|^2}$','$\mathsf{n_{total,\,low}}$',...
    '$\mathsf{n_{total,\,high}}$','location','northwest');
% l = legend('$\mathsf{\stackrel{-}{n}_{low}}$','$\mathsf{|\alpha_{0,\,low}|^2}$','$\mathsf{\stackrel{-}{n}_{high}}$',...
%     '$\mathsf{|\alpha_{0,\,high}|^2}$','$\mathsf{n_{total,\,low}}$',...
%     '$\mathsf{n_{total,\,high}}$','location','northwest');
l.Interpreter = 'latex';
l.FontSize = fontsize-2;
ylabel('Photon number');
xlabel([parameter ' (' xUnit ')']);
ylim([min([nCoherentLow nCoherentHigh])/10 max(meanNHigh)*10 ]);
if Pthr > 0
    xlim([0.1 20]);
    else 
    xlim([xmin max(Is)+10]);
end
graphicsSettings;
ax = gca;
ax.YMinorGrid = 'off';
set(ax,'FontSize',fontsize,'FontName',fontname);
savefig(['Husimiplots\PhotonNumbersFromHusimiFunctionsLowAndHigh' filenameOptions '.fig']);
print(['Husimiplots\PhotonNumbersFromHusimiFunctionsLowAndHigh' filenameOptions '.png'],'-dpng','-r300');
clf();

%% ratio of photon numbers 
if plotErrorbars
    errorbar(Is,nRatioLow,zeros(length(Is),1),nRatioErrLow,'ok','linewidth',1.5,'markerSize',msize);
    f = gca; f.XScale = 'log';f.YScale = 'log';hold on; 
    errorbar(Is,nRatioHigh,zeros(length(Is),1),nRatioErrHigh,'*k','linewidth',1.5,'markerSize',msize);
else
    semilogx(Is,nRatioLow,'ok',Is,nRatioHigh,'*k','markerSize',10,'markerSize',msize);
end
l = legend('Low','High','location','northwest');
l.FontSize = fontsize;
%ylabel('$|\alpha_{0}|^2/\bar n$','Interpreter','latex');
%ylabel('$\mathsf{|\alpha_{0}|^2/\stackrel{-}{n}}$','Interpreter','latex');
ylabel('$\mathsf{|\alpha_{0}|^2}/\bar \mathsf{n}$','Interpreter','latex','FontName','Arial');
xlabel([parameter ' (' xUnit ')']);
if Pthr > 0
    xlim([0.1 20]);
    else 
    xlim([xmin max(Is)+10]);
end
graphicsSettings;
ax = gca;
ax.YMinorGrid = 'off';
set(ax,'FontSize',fontsize,'FontName',fontname);
savefig(['Husimiplots\PhotonNumberRatioLowAndHigh' filenameOptions '.fig']);
print(['Husimiplots\PhotonNumberRatioLowAndHigh' filenameOptions '.png'],'-dpng','-r300');
clf();

%% g2
if plotErrorbars
    errorbar(Is,g2Low,g2ErrLow,'ok','linewidth',1.5,'markerSize',msize);
    f = gca; f.XScale = 'log';hold on; 
    errorbar(Is,g2High,g2ErrHigh,'*k','linewidth',1.5,'markerSize',msize);
else
    semilogx(Is,g2Low,'ok',Is,g2High,'*k','markerSize',msize);
end
l = legend('Low','High','location','southwest');
l.FontSize = fontsize;
ylabel('g^{(2)}(0)');
xlabel([parameter ' (' xUnit ')']);
ylim([1 2]);
if Pthr > 0
    xlim([0.1 20]);
     else 
    xlim([xmin max(Is)+10]);
end
graphicsSettings;
ax = gca;
set(ax,'FontSize',fontsize,'FontName',fontname);
savefig(['Husimiplots\g2LowAndHigh' filenameOptions '.fig']);
print(['Husimiplots\g2LowAndHigh' filenameOptions '.png'],'-dpng','-r300');
clf();

%% errors
plot(Is,nThermErr,'o',Is,nCohErr,'o',Is(CoherenceErr>0),CoherenceErr(CoherenceErr>0),'o');
legend('nThermErr','nCohErr','CoherenceErr');
savefig(['Errors' filenameOptions '.fig']);
print(['Errors' filenameOptions '.png'],'-dpng');

%% weights
semilogx(Is,weightLow,'o',Is,weightHigh,'o');
legend('Low','High','location','northwest');
ylabel('Relative abundance');
xlabel([parameter ' (' xUnit ')']);
graphicsSettings;
ax = gca;
set(ax,'FontSize',fontsize,'FontName',fontname);
savefig(['Husimiplots\PhotonNumbersFromHusimiFunctionsWeights' filenameOptions '.fig']);
print(['Husimiplots\PhotonNumbersFromHusimiFunctionsWeights' filenameOptions '.png'],'-dpng');
clf();

%% meanN

if plotErrorbars
    errorbar(Is,meanNLow,zeros(length(Is),1),meanNErrLow,'ok','linewidth',1.5,'markerSize',10);
    f = gca; f.XScale = 'log';f.YScale = 'log';hold on; 
    errorbar(Is,meanNHigh,zeros(length(Is),1),meanNErrHigh,'*k','linewidth',1.5,'markerSize',10);
else
   loglog(Is,meanNLow,'ok',Is,meanNHigh,'*k','markerSize',10,'markerSize',10);
end
l = legend('Low','High','location','northwest');
l.FontSize = fontsize;
ylabel('<n_{total}>');
xlabel([parameter ' (' xUnit ')']);
if Pthr > 0
    xlim([0.1 20]);
end
graphicsSettings;
ax = gca;
set(ax,'FontSize',fontsize,'FontName',fontname);
savefig(['Husimiplots\NTotalLowAndHigh' filenameOptions '.fig']);
print(['Husimiplots\NTotalLowAndHigh' filenameOptions '.png'],'-dpng','-r300');
clf();


end % function

