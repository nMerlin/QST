function [ dataStruct] = g22ChPhaseAveragedBinnedFromQuadratures( nResolution, channela, channelb, range, varargin)
%DLSERIES Batch processing of series of g2 measurements, with already
%computed quadratures that are in the folder 'Quadratures'
%Computes g2 either time resolved (set 'G2method', 'time') or photon number
%resolved (set 'G2method', 'bins').
%   'Weight': If set 'yes', g2 from the time resolved method is weighted 
%   with the respective photon numbers. 

% Optional input arguments
%% Validate and parse input arguments
p = inputParser;
defaultBins = 100; % to be used with the method 'bins'
addParameter(p,'Bins',defaultBins,@isnumeric);
parse(p,varargin{:});
c = struct2cell(p.Results);
[bins] = c{:};

%% Variables
dataStruct = struct('filename',{},'Delay', {},'NPhotonsCh1', {},'G2Ch1',{},'NPhotonsCh2',{},'G2Ch2',{},'G2Ch1Ch2',{});

%% Create data overview
dispstat('','init');
dispstat('Checking filenames ...','timestamp','keepthis');
Contents = dir('Quadratures');
name = {Contents.name};

cd('Quadratures');
for iStruct =  1:length(Contents) 
    %get filename
    filename = cell2mat(name(iStruct));
    if strcmp(filename,'.') || strcmp(filename,'..')
        continue
    end
    dataStruct(iStruct).filename = filename;
    % get current or power
    %currentToken = regexpi(filename,'-([0123456789.]*)mA','tokens');
%     currentToken = regexpi(filename,'mW-([0123456789.]*)mW','tokens');
     %currentToken = regexpi(filename,'([0123456789,]*)mW-5mW','tokens');
     %currentToken = regexpi(filename,'([0123456789,]*)mW','tokens'); 
     currentToken = regexpi(filename,'mat([0123456789,]*)','tokens');
     currentToken{1}=strrep(currentToken{1},',','.');
    dataStruct(iStruct).Delay = str2double(cell2mat(currentToken{1}));
    
    %%% Load data
    load(filename);
    
    % First channel
     switch channela
        case 1
            Xa = X1;
        case 2
            Xa = X2;
        case 3
            Xa = X3;
     end
    [~, ada,meang2, XOutA, IndicesA, nA] = g2Bins(Xa, nResolution, bins, 'xx','Plot','no');

    g2values = meang2;
    nPhotonsA = mean(ada(round(bins/3) : round(bins*2/3)));
    
    disp(['Delay: ',num2str(dataStruct(iStruct).Delay)]) 
    disp(['First g2 of Channel 1: ',num2str(g2values)])
    disp(['First nPhotons: ',num2str(nPhotonsA)])
    dataStruct(iStruct).G2Ch1 = g2values;
    dataStruct(iStruct).NPhotonsCh1 = nPhotonsA;
   

    % second channel
     switch channelb
        case 1
            Xb = X1;
        case 2
            Xb = X2;
        case 3
            Xb = X3;
    end
    [~, ada,meang2, XOutB, IndicesB, nB] = g2Bins(Xb, nResolution, bins, 'xx','Plot','no');

    g2values = meang2;
    nPhotonsB = mean(ada(round(bins/3) : round(bins*2/3)));
    
    disp(['Delay: ',num2str(dataStruct(iStruct).Delay)]) 
    disp(['First g2 of Channel 2: ',num2str(g2values)])
    disp(['First nPhotons: ',num2str(nPhotonsB)])
    dataStruct(iStruct).G2Ch2 = g2values;
    dataStruct(iStruct).NPhotonsCh2 = nPhotonsB;
    %save(strcat(num2str(i+1,'%02d'),'-LOwithSIG-DL-',num2str(Delays(i)),'-Ch2.mat'),'X2','piezoSign2','theta2','uniformX2','uniformTheta2');

    % correlation between both channels
    % chose only middle bins (nur ein Bin nehmen? Dann Mittelwerte?)
    Xa = Xa(:);
    Xb = Xb(:);
    Xab = Xa + Xb; 
    XOutA = Xa( nA>=(1-range)*nPhotonsA & nA<=(1+range)*nPhotonsA  &  nB>=(1-range)*nPhotonsB & nB<=(1+range)*nPhotonsB);
    XOutB = Xb( nA>=(1-range)*nPhotonsA & nA<=(1+range)*nPhotonsA  &  nB>=(1-range)*nPhotonsB & nB<=(1+range)*nPhotonsB );
    Xab = Xab( nA>=(1-range)*nPhotonsA & nA<=(1+range)*nPhotonsA  &  nB>=(1-range)*nPhotonsB & nB<=(1+range)*nPhotonsB );
    
%     XOutA = XOutA(:,round(bins/3) : round(bins*2/3));
%     IndicesA = IndicesA(:,round(bins/3) : round(bins*2/3));
%     XOutB = XOutB(:,round(bins/3) : round(bins*2/3));
%     IndicesB = IndicesB(:,round(bins/3) : round(bins*2/3));
%     XOutA = XOutA(:);
%     XOutA = XOutA(~isnan(XOutA));
%     IndicesA = IndicesA(:);
%     IndicesA = IndicesA(~isnan(XOutA));
%     XOutB = XOutB(:);
%     XOutB = XOutB(~isnan(XOutB));
%     IndicesB = IndicesB(:);
%     IndicesB = IndicesB(~isnan(XOutB));
%     [IndicesA,I]=sort(IndicesA);
%     XOutA = XOutA(I);
%     [IndicesB,I]=sort(IndicesB);
%     XOutB = XOutB(I);
%     minlength = min(length(XOutA),length(XOutB));
%     XOutA = XOutA(1:minlength);
%     XOutB = XOutB(1:minlength);
%     IndicesA = IndicesA(1:minlength);
%     IndicesB = IndicesB(1:minlength);
%     Xab = XOutA(IndicesA==IndicesB)+XOutB(IndicesA==IndicesB);
    g2values = g22Ch(XOutA, XOutB, Xab,length(XOutA),length(XOutB),length(Xab));
    disp(['Delay: ',num2str(dataStruct(iStruct).Delay)])
    disp(['First g2 for two channels: ',num2str(g2values)])
    dataStruct(iStruct).G2Ch1Ch2 = g2values;
    
end % iStruct
cd('..');
Delays = cell2mat({dataStruct.Delay});
NPhvsDLCh1 = cell2mat({dataStruct.NPhotonsCh1});
G2vsDLCh1 = cell2mat({dataStruct.G2Ch1});
NPhvsDLCh2 = cell2mat({dataStruct.NPhotonsCh2});
G2vsDLCh2 = cell2mat({dataStruct.G2Ch2});
G2vsDLCh1Ch2 = cell2mat({dataStruct.G2Ch1Ch2});

c = 299792458;
timedelay = (Delays-Delays(1))*2*10^-3 /c *10^12;% time in ps
save(['Binned-G2-Ch1Ch2-channel1-' num2str(channela) '-channel2-' num2str(channelb) '-range-' num2str(range) '.mat'],...
    'Delays','NPhvsDLCh1','G2vsDLCh1','NPhvsDLCh2','G2vsDLCh2','G2vsDLCh1Ch2','timedelay');
figure(1),plot(timedelay ,NPhvsDLCh1,timedelay,NPhvsDLCh2); 
figure(2),plot(timedelay,G2vsDLCh1,timedelayG2vsDLCh2);
figure(3),plot(timedelay,G2vsDLCh1Ch2);
fontsize =22;
graphicsSettings;
xlim([min(timedelay) max(timedelay)]);
xlabel('time delay $\tau$(ps)','FontSize',fontsize,'Interpreter','latex');
ylabel('$  g^{(2)}_{12}(\tau)  $','FontSize',fontsize, 'Interpreter','latex');
savefig(['g2Ch1Ch2vsDelay-range-' num2str(range) '.fig']);
print(['g2Ch1Ch2vsDelay-range-' num2str(range)],'-dpng');
end % function

