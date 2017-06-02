function [] = coherentSeries(varargin)
%This function evaluates a measurement series of coherent states with
%different photon numbers.
%
% 'errorbars': Calculate the mean of uncertainties over all measured
% piezo-segments and plot errorbars according to the standard deviation

% Optional input arguments
verbose = 0;
quiet = 'notquiet';
errorbars = 0;
if nargin > 0
    for i = 1:nargin
        eval([varargin{i} '=1;']);
    end
end
if verbose == 0
    quiet = 'quiet';
end

Norm = 1/sqrt(2);
%Norm ist the factor in the relation between the quadratures and the ladder
%operators: q = Norm*(a^{+} + a), p = Norm*i*(a^{+} - a)
%typical values are 1/sqrt(2) or 1/2. 

dataStruct = struct('filename',{},'powerLO', {},'meanN',{}, 'unc', {},...
    'delQ', {}, 'delP', {});
dataStructLOonly = struct('filename',{},'number',{});

%%% Create data overview
dispstat('','init',quiet);
dispstat('Checking filenames ...','timestamp','keepthis',quiet);
rawDataContents = dir('raw-data');

% LOwithLO-files
for name = {rawDataContents.name}
    % Loop only over *LOwithLO.raw files
    filename = cell2mat(name);
    if not(isempty(regexpi(filename,'LOwithLO.raw.','match')))...
            || isempty(regexpi(filename,'LOwithLO.raw','match'))
        continue
    end
    
    % Fetch number of measurement
    numberToken = regexpi(filename,'^([0123456789]*)-','tokens');
    number = str2double(cell2mat(numberToken{1}));
    dataStruct(number).filename = filename;
    
    % Get LO power
    powerToken = regexpi(filename,'-([0123456789.]*)mW','tokens');
    dataStruct(number).powerLO = str2double(cell2mat(powerToken{1}));
end

% LOonly-files
for name = {rawDataContents.name}
    % Loop only over *LOonly.raw files
    filename = cell2mat(name);
    if not(isempty(regexpi(filename,'LOonly.raw.','match')))...
            || isempty(regexpi(filename,'LOonly.raw','match'))
        continue
    end
    
    % Fetch number of measurement
    numberToken = regexpi(filename,'^([0123456789]*)-','tokens');
    number = str2double(cell2mat(numberToken{1}));
    dataStructLOonly(number).filename = filename;
    dataStructLOonly(number).number = number;
end
LOnumbers = cell2mat({dataStructLOonly.number});

for number = 1:size(dataStruct,2)

    filenameSIG = dataStruct(number).filename;
    if isempty(filenameSIG)
        continue
    end    
    
    % Use available quadrature dataset, if possible
    mainContents = dir();
    datasetExisting = 0;
    for name = {mainContents.name}
        mainName = cell2mat(name);
        if not(isempty(regexpi(mainName,['quadratureDataset-' strrep(num2str(filenameSIG),...
                'LOwithLO.raw','')],'match')))
            datasetExisting = 1;
            M = load(mainName);
            X = M.X;
            theta = M.theta;
        else
            continue
        end
    end
    
    if datasetExisting == 0
        %find adequate LOonly-file
        LOnumber = max(LOnumbers(LOnumbers<=number));    
        filenameLO = dataStructLOonly(LOnumber).filename;

        dispstat(['mainPrepareData number ' num2str(number)],...
            'timestamp','keepthis',quiet);
        [ X, theta ] = mainPrepareData( filenameLO, filenameSIG );
    end
    
    dispstat(['discretize data ' num2str(number)],'timestamp','keepthis',quiet);
    [ Xdis, thetadis ] = discretizeTheta( X, theta, 100 );
    dispstat(['compute and plot Expectations ' num2str(number)],'timestamp','keepthis',quiet);
    [ ~,~, ~, ~, delQ, delP, meanUnc, ~, meanN, ~ ] =...
        computeExpectations2( Xdis, thetadis, strrep(num2str(filenameSIG),'-LOwithLO.raw','') );
    dataStruct(number).meanN = meanN;
    dataStruct(number).unc = meanUnc;
    dataStruct(number).delQ = delQ;
    dataStruct(number).delP = delP;   
end

% Modify dataStruct with NaN-Values to convert it to arrays
nPiezoSegments = 0;
for iStruct = 1:length(dataStruct)
    if length(dataStruct(iStruct).meanN) > nPiezoSegments
        nPiezoSegments = length(dataStruct(iStruct).meanN);
    end
end
for iStruct = 1:length(dataStruct)
    if isempty(dataStruct(iStruct).meanN)
        continue
    end
    meanN = dataStruct(iStruct).meanN;
    dataStruct(iStruct).meanN = ...
        [meanN; NaN(nPiezoSegments-length(meanN),1)];
    unc = dataStruct(iStruct).unc;
    dataStruct(iStruct).unc = ...
        [unc; NaN(nPiezoSegments-length(unc),1)];
    delQ = dataStruct(iStruct).delQ;
    dataStruct(iStruct).delQ = ...
        [delQ; NaN(nPiezoSegments-length(delQ),1)];
    delP = dataStruct(iStruct).delP;
    dataStruct(iStruct).delP = ...
        [delP; NaN(1,nPiezoSegments-length(delP),1)];
end
meanNs = cell2mat({dataStruct.meanN});
uncs = cell2mat({dataStruct.unc});
delQs = cell2mat({dataStruct.delQ});
delPs = cell2mat({dataStruct.delP});


dispstat('Plot uncertainties over mean photon number','timestamp',...
    'keepthis',quiet);
if errorbars == 0 % Create a separate plot for each piezo segment
    for iSegment = 1:nPiezoSegments
        close all;
        semilogx(meanNs(iSegment,:),delQs(iSegment,:),'o',meanNs(iSegment,:),delPs(iSegment,:),'o', meanNs(iSegment,:),uncs(iSegment,:), 'o');
        hold on;
        semilogx(1:max(ceil(meanNs(iSegment,:))),Norm^2*ones(1,max(ceil(meanNs(iSegment,:)))),'k-');
        legend( '$\Delta Q$','$\Delta P$', '$\Delta Q \cdot \Delta P$ ','Location','northwest');
        axis ([1 max(meanNs(iSegment,:))+1 0 1.5]);
        set(0,'DefaultLegendInterpreter','latex');
        set(0,'DefaultTextInterpreter','latex');
        xlabel('mean Photonnumber');
        ylabel('Uncertainty [a. u.]');    
        print(strcat('UncertaintiesOverPhotonnumber-Segment',num2str(iSegment)),'-dpng');
        hold off;
    end
else % Create one plot that accounts for all piezo segments with error bars
    close all;
    hold on;
    %stdMeanNs = nanstd(meanNs);
    meanNs = nanmean(meanNs);
    stdDelQs = nanstd(delQs);
    delQs = nanmean(delQs);
    stdDelPs = nanstd(delPs);
    delPs = nanmean(delPs);
    stdUncs = nanstd(uncs);
    uncs = nanmean(uncs);
    errorbar(meanNs,delQs,stdDelQs,'o');
    errorbar(meanNs,delPs,stdDelPs,'o');
    errorbar(meanNs,uncs,stdUncs,'o');
    % Horizontal error bars are usually smaller than the dot
    % herrorbar(meanNs,uncs,stdMeanNs,stdMeanNs,'o');
    set(gca, 'xscale', 'log');
    semilogx(1:max(ceil(meanNs)),Norm^2*ones(1,max(ceil(meanNs))),'k-');
    legend( '$\Delta Q$','$\Delta P$', '$\Delta Q \cdot \Delta P$ ','Location','northwest');
    axis ([1 max(ceil(meanNs)) 0 1.5]);
    set(0,'DefaultLegendInterpreter','latex');
    set(0,'DefaultTextInterpreter','latex');
    xlabel('mean Photonnumber');
    ylabel('Uncertainty [a. u.]');
    print('UncertaintiesOverPhotonnumber-ErrorBars','-dpng');
    hold off;
end

% Saving important results
dateString = datestr(datetime('now'),'yyyymmddTHHMMSS');
save(['photonnumberDataset-' dateString '.mat'], ...
    'meanNs', 'delQs', 'delPs', 'uncs');

end