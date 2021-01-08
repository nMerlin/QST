function HusimiSeriesFromQuadratures( chAssign, varargin)
%chAssign... 

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
defaultFitMethod = 'NLSQ-LAR';
addParameter(p,'FitMethod',defaultFitMethod,@isstr);
defaultScale = true; %scales the O1 and O2 so they have the same photon number
addParameter(p,'Scale',defaultScale,@islogical);
parse(p,varargin{:});
c = struct2cell(p.Results);
[fitMethod,loadExistent,LOpower,parameter,plotErrorbars,range,saveHusimi,scale,showLegend,xUnit] = c{:};

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
            theta = reshape(XpsFast,[a*b c]); %not real theta, but needed 
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
        [H, binsO1, binsO2,r,nTherm,nThermErr,nCoherent,nCohErr,meanN,Coherence,CoherenceErr] = ...
            plotHusimiAndCut(O1,O2,0,0,0,'Filename',['Husimiplots\' filename ...
            '-scaled-' num2str(scaled) '-showLegend-' num2str(showLegend)],...
            'ShowLegend',showLegend,'FitMethod',fitMethod);
         close all; 
        if saveHusimi
            save(['mat-data\' filename],'H','binsO1','binsO2','-append');
        end
        rscMax = r;
        [nThermMax,nThermHigh,nThermLow] = deal(nTherm);
        [nCoherentMax,nCoherentHigh,nCoherentLow] = deal(nCoherent);
        [meanNMax,meanNHigh,meanNLow] = deal(meanN);
        [weightLow,weightHigh] = deal(0);
        [CoherenceMax,CoherenceHigh,CoherenceLow] = deal(Coherence);
        [nThermErrMax,nThermErrHigh,nThermErrLow] = deal(nThermErr);
        [nCohErrMax,nCohErrHigh,nCohErrLow] = deal(nCohErr);
        [CoherenceErrMax,CoherenceErrHigh,CoherenceErrLow] = deal(CoherenceErr);
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
        [~, ~, ~,rLow,nThermLow,nThermErrLow,nCoherentLow,nCohErrLow,meanNLow,CoherenceLow,CoherenceErrLow] = ...
            plotHusimiAndCut(O1s,O2s,0,0,0,'Filename',...
            ['Husimiplots\' filename '-scaled-' num2str(scaled) '-showLegend-' num2str(showLegend) '-LowPhotonNumber'],...
            'ShowLegend',showLegend,'FitMethod',fitMethod);
         close all;
         % for the high one 
        iSelHigh = find(nPsFastVec >= min(rangeHigh) & nPsFastVec <= max(rangeHigh));
        O1s = O1(iSelHigh);
        O2s = O2(iSelHigh);
        [~, ~, ~,rHigh,nThermHigh,nThermErrHigh,nCoherentHigh,nCohErrHigh,meanNHigh,CoherenceHigh,CoherenceErrHigh] = ...
            plotHusimiAndCut(O1s,O2s,0,0,0,'Filename',...
            ['Husimiplots\' filename '-scaled-' num2str(scaled) '-showLegend-' num2str(showLegend) '-HighPhotonNumber'],...
            'ShowLegend',showLegend,'FitMethod',fitMethod);
         close all;
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

save(['Husimiplots\HusimiResultsScaled-FitMethod-' fitMethod '-scaled-' num2str(scaled) '.mat'],'Is','r','nTherm','meanN','nCoherent',...
    'rMax','nThermMax','meanNMax','nCoherentMax','weightLow', 'weightHigh', 'nThermHigh', 'nCoherentHigh', 'meanNHigh', 'nThermLow', 'meanNLow','nCoherentLow',...
    'CoherenceHigh','CoherenceMax','Coherence','CoherenceLow','nThermErr','nCohErr','CoherenceErr','nThermErrMax','nCohErrMax','CoherenceErrMax');
xlswrite(['Husimiplots\HusimiResultsScaled-' fitMethod '-scaled-' num2str(scaled) '.xls'],[Is' r' nTherm' meanN' nCoherent' ...
    rMax' nThermMax' meanNMax' nCoherentMax' weightLow' weightHigh' nThermHigh' nCoherentHigh' meanNHigh' nThermLow' meanNLow' nCoherentLow' ...
     CoherenceHigh' CoherenceMax' Coherence' CoherenceLow' nThermErr' nCohErr' CoherenceErr' nThermErrMax' nCohErrMax' CoherenceErrMax']);

 %% plotting
if strcmp(parameter,'Power')
    parameter = 'Excitation power P';
end
fontsize = 20;
fontname = 'Arial';
filenameOptions = ['-scaled-' num2str(scale) '-FitMethod-' fitMethod '-errorbars-' num2str(plotErrorbars)]; 

%% coherence
if plotErrorbars
    errorbar(Is,CoherenceLow,zeros(length(Is),1),CoherenceErr,'ok','linewidth',1.5,'markerSize',10);
    f = gca;f.XScale = 'log'; hold on; 
    errorbar(Is,CoherenceHigh,zeros(length(Is),1),CoherenceErr,'*k','linewidth',1.5,'markerSize',10);
else
    semilogx(Is,CoherenceLow,'ok',Is,CoherenceHigh,'*k','markerSize',10);
end
l=legend('$\mathcal{C}_\mathrm{low}$','$\mathcal{C}_\mathrm{high}$','location',...
    'northwest');
l.Interpreter = 'latex';
l.FontSize = fontsize+8;
ylabel('Coherence');
xlabel([parameter ' (' xUnit ')']);
graphicsSettings;
ax = gca;
set(ax,'FontSize',fontsize,'FontName',fontname);
savefig(['Husimiplots\CoherenceFromHusimiFunctions-semilogx' filenameOptions '.fig']);
print(['Husimiplots\CoherenceFromHusimiFunctions-semilogx' filenameOptions '.png'],'-dpng','-r700');
clf();

%% photon numbers all together
if plotErrorbars
    errorbar(Is,nThermLow,zeros(length(Is),1),nThermErr,'ok','linewidth',1.5);
    f = gca; f.XScale = 'log';f.YScale = 'log';hold on; 
    errorbar(Is,nCoherentLow,zeros(length(Is),1),nCohErr,'or','linewidth',1.5);
    errorbar(Is,nThermHigh,zeros(length(Is),1),nThermErr,'*k','linewidth',1.5);
    errorbar(Is,nCoherentHigh,zeros(length(Is),1),nCohErr,'*r','linewidth',1.5);
else
    loglog(Is,nThermLow,'ok',Is,nCoherentLow,'or',Is,nThermHigh,'*k',Is,nCoherentHigh,'*r','markerSize',7);
end
l = legend('$\bar n_\mathrm{low}$','$|\alpha_{0,\mathrm{low}}|^2$','$\bar n_\mathrm{high}$',...
    '$|\alpha_{0,\mathrm{high}}|^2$','location','northwest');
l.Interpreter = 'latex';
l.FontSize = fontsize+8;
ylabel('Photon number');
xlabel([parameter ' (' xUnit ')']);
graphicsSettings;
ax = gca;
set(ax,'FontSize',fontsize,'FontName',fontname);
savefig(['Husimiplots\PhotonNumbersFromHusimiFunctionsLowAndHigh' filenameOptions '.fig']);
print(['Husimiplots\PhotonNumbersFromHusimiFunctionsLowAndHigh' filenameOptions '.png'],'-dpng','-r700');
clf();

%% errors
plot(Is,nThermErr,'o',Is,nCohErr,'o',Is(CoherenceErr>0),CoherenceErr(CoherenceErr>0),'o');
legend('nThermErr','nCohErr','CoherenceErr');
savefig(['Errors' filenameOptions '.fig']);
print(['Errors' filenameOptions '.png'],'-dpng');

%%
if plotErrorbars
    errorbar(Is,nTherm,nThermErr,'o','linewidth',1.5);
    f = gca; f.XScale = 'log';f.YScale = 'log';hold on; 
    errorbar(Is,nCoherent,nCohErr,'o','linewidth',1.5);
else
    loglog(Is,nTherm,'o');hold on;loglog(Is,nCoherent,'o');
end
legend('n_{Thermal}','n_{Coherent}','location','northwest');
ylabel('Photon number');
xlabel([parameter ' (' xUnit ')']);
graphicsSettings;
ax = gca;
set(ax,'FontSize',fontsize,'FontName',fontname);
savefig(['Husimiplots\PhotonNumbersFromHusimiFunctionsWeighted' filenameOptions '.fig']);
print(['Husimiplots\PhotonNumbersFromHusimiFunctionsWeighted' filenameOptions '.png'],'-dpng');
clf();

loglog(Is,nCoherent./nTherm,'o');
legend('n_{Coherent}/n_{Thermal}','location','northwest');
ylabel('Ratio');
xlabel([parameter ' (' xUnit ')']);
graphicsSettings;
ax = gca;
set(ax,'FontSize',fontsize,'FontName',fontname);
savefig(['Husimiplots\PhotonNumberRatioFromHusimiFunctionsWeighted' filenameOptions '.fig']);
print(['Husimiplots\PhotonNumberRatioFromHusimiFunctionsWeighted' filenameOptions '.png'],'-dpng');
clf();

semilogx(Is,nCoherent./nTherm,'o');
legend('n_{Coherent}/n_{Thermal}','location','northwest');
ylabel('Ratio');
xlabel([parameter ' (' xUnit ')']);
graphicsSettings;
ax = gca;
set(ax,'FontSize',fontsize,'FontName',fontname);
savefig(['Husimiplots\PhotonNumberRatioFromHusimiFunctionsWeighted-semilogx' filenameOptions '.fig']);
print(['Husimiplots\PhotonNumberRatioFromHusimiFunctionsWeighted-semilogx' filenameOptions '.png'],'-dpng');
clf();

loglog(Is,nThermMax,'o');hold on;loglog(Is,nCoherentMax,'o');
legend('n_{Thermal}','n_{Coherent}','location','northwest');
ylabel('Photon number');
xlabel([parameter ' (' xUnit ')']);
graphicsSettings;
ax = gca;
set(ax,'FontSize',fontsize,'FontName',fontname);
savefig(['Husimiplots\PhotonNumbersFromHusimiFunctionsMax-FitMethod' filenameOptions '.fig']);
print(['Husimiplots\PhotonNumbersFromHusimiFunctionsMax-FitMethod' filenameOptions '.png'],'-dpng');
clf();

loglog(Is,nCoherentMax./nThermMax,'o');
legend('n_{Coherent}/n_{Thermal}','location','northwest');
ylabel('Ratio');
xlabel([parameter ' (' xUnit ')']);
graphicsSettings;
ax = gca;
set(ax,'FontSize',fontsize,'FontName',fontname);
savefig(['Husimiplots\PhotonNumberRatioFromHusimiFunctionsMax-FitMethod' filenameOptions '.fig']);
print(['Husimiplots\PhotonNumberRatioFromHusimiFunctionsMax-FitMethod' filenameOptions '.png'],'-dpng');
clf();

semilogx(Is,nThermLow,'o',Is,nCoherentLow,'o',Is,nThermHigh,'o',Is,nCoherentHigh,'o');
legend('n_{Thermal,Low}','n_{Coherent,Low}','n_{Thermal,High}','n_{Coherent,High}','location','northwest');
ylabel('Photon number');
xlabel([parameter ' (' xUnit ')']);
graphicsSettings;
ax = gca;
set(ax,'FontSize',fontsize,'FontName',fontname);
savefig(['Husimiplots\PhotonNumbersFromHusimiFunctionsLowAndHigh-semilogx' filenameOptions '.fig']);
print(['Husimiplots\PhotonNumbersFromHusimiFunctionsLowAndHigh-semilogx' filenameOptions '.png'],'-dpng');
clf();

loglog(Is,meanNLow,'o',Is,meanNHigh,'o');
legend('n_{total,Low}','n_{total,High}','location','northwest');
ylabel('Photon number');
xlabel([parameter ' (' xUnit ')']);
graphicsSettings;
ax = gca;
set(ax,'FontSize',fontsize,'FontName',fontname);
savefig(['Husimiplots\PhotonNumbersFromHusimiFunctionsMeanNLowAndHigh' filenameOptions '.fig']);
print(['Husimiplots\PhotonNumbersFromHusimiFunctionsMeanNLowAndHigh' filenameOptions '.png'],'-dpng');
clf();

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


end % function

