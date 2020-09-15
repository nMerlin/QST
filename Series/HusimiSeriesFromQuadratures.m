function HusimiSeriesFromQuadratures( chAssign, varargin)
%chAssign... 

% Optional input arguments
%% Validate and parse input arguments
p = inputParser;
defaultParameter = 'power'; % which parameter was changed during the series
addParameter(p,'Parameter',defaultParameter);
defaultXUnit = 'mW'; % unit  of the parameter
addParameter(p,'XUnit',defaultXUnit);
parse(p,varargin{:});
c = struct2cell(p.Results);
[parameter,xUnit] = c{:};

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
    
    
%     try
%         load(['mat-data\' filename],'Hsc','binsO1sc','binsO2sc');
%          [Hsc, binsO1sc, binsO2sc,rsc,nThermsc,nCoherentsc,meanNsc] = plotHusimiAndCut(0,0,Hsc,binsO1sc,binsO2sc,'Filename',['Husimiplots\scaled-' filename]);
%          close all;
%         save(['mat-data\' filename],'Hsc','binsO1sc','binsO2sc','-append');
%     catch
        load(['mat-data\' filename]);
        quadratures = zeros([size(X1) 3]);
        quadratures(:,:,:,1) = X1;
        quadratures(:,:,:,2) = X2;
        quadratures(:,:,:,3) = X3;
        XpsFast = quadratures(:,:,:,chAssign(2));
        XpsSlow = quadratures(:,:,:,chAssign(3));
        clear('quadratures');
        [~,nPsFast,nPsSlow] = nPhotons(XpsFast,XpsFast,XpsSlow);
        
        if ~exist('O1','var')  % make orthogonal quadratures 
            [a,b,c] = size(XpsFast);
            theta = reshape(XpsFast,[a*b c]); %not real theta, but needed 
            [O1,O2,~,~,iOrth] = selectOrthogonal(XpsFast,XpsSlow,XpsSlow,theta,piezoSign,'width',0.05);
            save(['mat-data\' filename],'O1','O2','iOrth','-append');
        end

        [H, binsO1, binsO2] = plotHusimiAndCut(O1,O2,0,0,0,'Filename',['Husimiplots\' filename],'MakeFit',false );
        close all;

        O1scaled = O1*sqrt(mean([nPsFast nPsSlow])/nPsFast);
        O2scaled = O2*sqrt(mean([nPsFast nPsSlow])/nPsSlow);
        [Hsc, binsO1sc, binsO2sc,rsc,nThermsc,nCoherentsc,meanNsc] = ...
            plotHusimiAndCut(O1scaled,O2scaled,0,0,0,'Filename',['Husimiplots\scaled-' filename],'MakeFit',true);
         close all;        
        save(['mat-data\' filename],'H','Hsc','binsO1','binsO1sc','binsO2','binsO2sc','-append');
        clear X1 X2 X3 theta piezoSign O1 O2 'H' 'Hsc' 'binsO1' 'binsO1sc' 'binsO2' 'binsO2sc' O1scaled O2scaled 
        clear XpsFast XpsSlow nPsFast nPsSlow
%     end
    
              
    dataStruct(iStruct).r = rsc;
    dataStruct(iStruct).nTherm = nThermsc;
    dataStruct(iStruct).nCoherent = nCoherentsc;
    dataStruct(iStruct).meanN = meanNsc;
    
end % iStruct

Is = cell2mat({dataStruct.I});
r = cell2mat({dataStruct.r});
nTherm = cell2mat({dataStruct.nTherm});
nCoherent = cell2mat({dataStruct.nCoherent});
meanN = cell2mat({dataStruct.meanN});

save('Husimiplots\HusimiResultsScaled.mat','Is','r','nTherm','meanN','nCoherent');
xlswrite('Husimiplots\HusimiResultsScaled.xls',[Is' r' nTherm' meanN' nCoherent']);

loglog(Is,nTherm,'o');hold on;loglog(Is,nCoherent,'o');loglog(Is,meanN,'o');
legend('n_{Thermal}','n_{Coherent}','n_{total}','location','best');
ylabel('photon number');
xlabel([parameter ' (' xUnit ')']);
graphicsSettings;
savefig('Husimiplots\PhotonNumbersFromHusimiFunctions.fig');
print('Husimiplots\PhotonNumbersFromHusimiFunctions.png','-dpng');

end % function

