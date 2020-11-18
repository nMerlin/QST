function HusimiSeriesFromQuadratures( chAssign, varargin)
%chAssign... 

% Optional input arguments
%% Validate and parse input arguments
p = inputParser;
defaultParameter = 'power'; % which parameter was changed during the series
addParameter(p,'Parameter',defaultParameter);
defaultXUnit = 'mW'; % unit  of the parameter
addParameter(p,'XUnit',defaultXUnit);
defaultRange = 0.5; % if maxN-minN > range*medN, there are two fits made 
%for higher and for lower photon numbers.  
addParameter(p,'range',defaultRange);
parse(p,varargin{:});
c = struct2cell(p.Results);
[parameter,range,xUnit] = c{:};

%% Variables
dataStruct = struct('filename',{},'I',{},'r',{},'nTherm',{},'nCoherent',{},'meanN',{});

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
        case 'power'    
            %currentToken = regexpi(filename,'([0123456789,]*)mW','tokens');
            currentToken = regexpi(filename,'([0123456789,]*)mW-4mW','tokens');
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
    dispstat(['load ' filename],...
        'timestamp','keepthis','notquiet');   
    
        load(['mat-data\' filename]);
        
        if ~exist('nPsFast','var')     
            load(['mat-data\' filename]);
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
  
        
        maxN = max(nPsFastVec(:));
        minN = min(nPsFastVec(:));
        medN = mean([maxN minN]);
        if maxN-minN < range*medN  || nPsFast <1     
        % If the fluctuation is small, there is probably only one state all the
        % time, so one fit is enough. 
            [H, binsO1, binsO2] = plotHusimiAndCut(O1,O2,0,0,0,'Filename',['Husimiplots\' filename],'MakeFit',false );
            close all;
            O1scaled = O1*sqrt(mean([nPsFast nPsSlow])/nPsFast);
            O2scaled = O2*sqrt(mean([nPsFast nPsSlow])/nPsSlow);
            [Hsc, binsO1sc, binsO2sc,rsc,nThermsc,nThermErr,nCoherentsc,nCohErr,meanNsc,Coherence] = ...
                plotHusimiAndCut(O1scaled,O2scaled,0,0,0,'Filename',['Husimiplots\scaled-' filename],'MakeFit',true);
             close all;        
            save(['mat-data\' filename],'H','Hsc','binsO1','binsO1sc','binsO2','binsO2sc','-append');
            rscMax = rsc;
            [nThermscMax,nThermscHigh,nThermscLow] = deal(nThermsc);
            [nCoherentscMax,nCoherentscHigh,nCoherentscLow] = deal(nCoherentsc);
            [meanNscMax,meanNscHigh,meanNscLow] = deal(meanNsc);
            [weightLow,weightHigh] = deal(0);
            [CoherenceMax,CoherenceHigh,CoherenceLow] = deal(Coherence);
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
            [H, binsO1, binsO2] = plotHusimiAndCut(O1s,O2s,0,0,0,'Filename',['Husimiplots\' filename '-LowPhotonNumber'],'MakeFit',false );
            close all;
            O1scaled = O1s*sqrt(mean([nPsFast nPsSlow])/nPsFast);
            O2scaled = O2s*sqrt(mean([nPsFast nPsSlow])/nPsSlow);
            [Hsc, binsO1sc, binsO2sc,rscLow,nThermscLow,nThermErrLow,nCoherentscLow,nCohErrLow,meanNscLow,CoherenceLow] = ...
                plotHusimiAndCut(O1scaled,O2scaled,0,0,0,'Filename',['Husimiplots\scaled-' filename '-LowPhotonNumber'],'MakeFit',true);
             close all;
             % for the high one 
            iSelHigh = find(nPsFastVec >= min(rangeHigh) & nPsFastVec <= max(rangeHigh));
            O1s = O1(iSelHigh);
            O2s = O2(iSelHigh);
            [H, binsO1, binsO2] = plotHusimiAndCut(O1s,O2s,0,0,0,'Filename',['Husimiplots\' filename '-HighPhotonNumber'],'MakeFit',false );
            close all;
            O1scaled = O1s*sqrt(mean([nPsFast nPsSlow])/nPsFast);
            O2scaled = O2s*sqrt(mean([nPsFast nPsSlow])/nPsSlow);
            [Hsc, binsO1sc, binsO2sc,rscHigh,nThermscHigh,nThermErrHigh,nCoherentscHigh,nCohErrHigh,meanNscHigh,CoherenceHigh] = ...
                plotHusimiAndCut(O1scaled,O2scaled,0,0,0,'Filename',['Husimiplots\scaled-' filename '-HighPhotonNumber'],'MakeFit',true);
             close all;
            % their respective weights
            weightLow = length(iSelLow)/(length(iSelLow)+length(iSelHigh));
            weightHigh = length(iSelHigh)/(length(iSelLow)+length(iSelHigh));
            % weighted results
            rsc = weightLow*rscLow + weightHigh*rscHigh;
            nThermsc = weightLow*nThermscLow + weightHigh*nThermscHigh;
            nCoherentsc = weightLow*nCoherentscLow + weightHigh*nCoherentscHigh;
            meanNsc = weightLow*meanNscLow + weightHigh*meanNscHigh;
            Coherence = weightLow*CoherenceLow + weightHigh*CoherenceHigh;
            % maximum results
            rscMax = max([rscLow rscHigh]);
            nThermscMax = max([nThermscLow nThermscHigh]);
            nCoherentscMax = max([nCoherentscLow nCoherentscHigh]);
            meanNscMax = max([meanNscLow meanNscHigh]);
            CoherenceMax = max([CoherenceLow CoherenceHigh]);
        end            
%     end  
    clear X1 X2 X3 theta piezoSign O1 O2 'H' 'Hsc' 'binsO1' 'binsO1sc' 'binsO2' 'binsO2sc' O1scaled O2scaled 
    clear XpsFast XpsSlow nPsFast nPsSlow nPsFastVec
    dataStruct(iStruct).r = rsc;
    dataStruct(iStruct).nTherm = nThermsc;
    dataStruct(iStruct).nCoherent = nCoherentsc;
    dataStruct(iStruct).meanN = meanNsc;  
    dataStruct(iStruct).Coherence = Coherence;  
    dataStruct(iStruct).rMax = rscMax;
    dataStruct(iStruct).nThermMax = nThermscMax;
    dataStruct(iStruct).nCoherentMax = nCoherentscMax;
    dataStruct(iStruct).meanNMax = meanNscMax; 
    dataStruct(iStruct).CoherenceMax = CoherenceMax;  
    dataStruct(iStruct).weightLow = weightLow;
    dataStruct(iStruct).weightHigh = weightHigh;
    dataStruct(iStruct).nThermHigh = nThermscHigh;
    dataStruct(iStruct).nCoherentHigh = nCoherentscHigh;
    dataStruct(iStruct).meanNHigh = meanNscHigh; 
    dataStruct(iStruct).CoherenceHigh = CoherenceHigh;  
    dataStruct(iStruct).nThermLow = nThermscLow;
    dataStruct(iStruct).nCoherentLow = nCoherentscLow;
    dataStruct(iStruct).meanNLow = meanNscLow;
    dataStruct(iStruct).CoherenceLow = CoherenceLow;  
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

save('Husimiplots\HusimiResultsScaled.mat','Is','r','nTherm','meanN','nCoherent',...
    'rMax','nThermMax','meanNMax','nCoherentMax','weightLow', 'weightHigh', 'nThermHigh', 'nCoherentHigh', 'meanNHigh', 'nThermLow', 'meanNLow','nCoherentLow',...
    'CoherenceHigh','CoherenceMax','Coherence','CoherenceLow');
xlswrite('Husimiplots\HusimiResultsScaled.xls',[Is' r' nTherm' meanN' nCoherent' ...
    rMax' nThermMax' meanNMax' nCoherentMax' weightLow' weightHigh' nThermHigh' nCoherentHigh' meanNHigh' nThermLow' meanNLow' nCoherentLow' ...
     CoherenceHigh' CoherenceMax' Coherence' CoherenceLow']);

 
semilogx(Is,CoherenceLow,'o',Is,CoherenceHigh,'o',Is,Coherence,'o');
legend('Coherence_{Low}','Coherence_{High}','Coherence_{Weighted}','location','northwest');
ylabel('Coherence');
xlabel([parameter ' (' xUnit ')']);
graphicsSettings;
savefig('Husimiplots\CoherenceFromHusimiFunctions-semilogx.fig');
print('Husimiplots\CoherenceFromHusimiFunctions-semilogx.png','-dpng');
clf();
 
loglog(Is,nTherm,'o');hold on;loglog(Is,nCoherent,'o');
legend('n_{Thermal}','n_{Coherent}','location','northwest');
ylabel('photon number');
xlabel([parameter ' (' xUnit ')']);
graphicsSettings;
savefig('Husimiplots\PhotonNumbersFromHusimiFunctionsWeighted.fig');
print('Husimiplots\PhotonNumbersFromHusimiFunctionsWeighted.png','-dpng');
clf();

loglog(Is,nCoherent./nTherm,'o');
legend('n_{Coherent}/n_{Thermal}','location','northwest');
ylabel('ratio');
xlabel([parameter ' (' xUnit ')']);
graphicsSettings;
savefig('Husimiplots\PhotonNumberRatioFromHusimiFunctionsWeighted.fig');
print('Husimiplots\PhotonNumberRatioFromHusimiFunctionsWeighted.png','-dpng');
clf();

semilogx(Is,nCoherent./nTherm,'o');
legend('n_{Coherent}/n_{Thermal}','location','northwest');
ylabel('ratio');
xlabel([parameter ' (' xUnit ')']);
graphicsSettings;
savefig('Husimiplots\PhotonNumberRatioFromHusimiFunctionsWeighted-semilogx.fig');
print('Husimiplots\PhotonNumberRatioFromHusimiFunctionsWeighted-semilogx.png','-dpng');
clf();

loglog(Is,nThermMax,'o');hold on;loglog(Is,nCoherentMax,'o');
legend('n_{Thermal}','n_{Coherent}','location','northwest');
ylabel('photon number');
xlabel([parameter ' (' xUnit ')']);
graphicsSettings;
savefig('Husimiplots\PhotonNumbersFromHusimiFunctionsMax.fig');
print('Husimiplots\PhotonNumbersFromHusimiFunctionsMax.png','-dpng');
clf();

loglog(Is,nCoherentMax./nThermMax,'o');
legend('n_{Coherent}/n_{Thermal}','location','northwest');
ylabel('ratio');
xlabel([parameter ' (' xUnit ')']);
graphicsSettings;
savefig('Husimiplots\PhotonNumberRatioFromHusimiFunctionsMax.fig');
print('Husimiplots\PhotonNumberRatioFromHusimiFunctionsMax.png','-dpng');
clf();

loglog(Is,nThermLow,'o',Is,nCoherentLow,'o',Is,nThermHigh,'o',Is,nCoherentHigh,'o');
legend('n_{Thermal,Low}','n_{Coherent,Low}','n_{Thermal,High}','n_{Coherent,High}','location','northwest');
ylabel('photon number');
xlabel([parameter ' (' xUnit ')']);
graphicsSettings;
savefig('Husimiplots\PhotonNumbersFromHusimiFunctionsLowAndHigh.fig');
print('Husimiplots\PhotonNumbersFromHusimiFunctionsLowAndHigh.png','-dpng');
clf();

semilogx(Is,nThermLow,'o',Is,nCoherentLow,'o',Is,nThermHigh,'o',Is,nCoherentHigh,'o');
legend('n_{Thermal,Low}','n_{Coherent,Low}','n_{Thermal,High}','n_{Coherent,High}','location','northwest');
ylabel('photon number');
xlabel([parameter ' (' xUnit ')']);
graphicsSettings;
savefig('Husimiplots\PhotonNumbersFromHusimiFunctionsLowAndHigh-semilogx.fig');
print('Husimiplots\PhotonNumbersFromHusimiFunctionsLowAndHigh-semilogx.png','-dpng');
clf();

loglog(Is,meanNLow,'o',Is,meanNHigh,'o');
legend('n_{total,Low}','n_{total,High}','location','northwest');
ylabel('photon number');
xlabel([parameter ' (' xUnit ')']);
graphicsSettings;
savefig('Husimiplots\PhotonNumbersFromHusimiFunctionsMeanNLowAndHigh.fig');
print('Husimiplots\PhotonNumbersFromHusimiFunctionsMeanNLowAndHigh.png','-dpng');
clf();

semilogx(Is,weightLow,'o',Is,weightHigh,'o');
legend('Low','High','location','northwest');
ylabel('relative abundance');
xlabel([parameter ' (' xUnit ')']);
graphicsSettings;
savefig('Husimiplots\PhotonNumbersFromHusimiFunctionsWeights.fig');
print('Husimiplots\PhotonNumbersFromHusimiFunctionsWeights.png','-dpng');
clf();


end % function

