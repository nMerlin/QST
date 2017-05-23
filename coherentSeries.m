function [] = coherentSeries()
%This function evaluates a measurement series of coherent states with
%different photon numbers.
% Optional input arguments
verbose = 0;
quiet = 'notquiet';
if nargin > 0
    for i = 1:nargin
        eval([varargin{i} '=1;']);
    end
end
if verbose == 0
    quiet = 'quiet';
end

dataStruct = struct('filename',{},'powerLO', {},'meanN',{}, 'unc', {}, 'DelQ', {}, 'DelP', {});
dataStructLOonly = struct('filename',{},'number',{});

%%% Create data overview
dispstat('','init',quiet);
dispstat('Checking filenames ...','timestamp','keepthis',quiet);
rawDataContents = dir('raw-data');

%LOwithLO-files
for name = {rawDataContents.name}
    % Loop only over *LOwithLO.raw files
    filename = cell2mat(name);
    if not(isempty(regexpi(filename,'LOwithLO.raw.','match'))) || isempty(regexpi(filename,'LOwithLO.raw','match'))
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

%LOonly-files
for name = {rawDataContents.name}
    % Loop only over *LOonly.raw files
    filename = cell2mat(name);
    if not(isempty(regexpi(filename,'LOonly.raw.','match'))) || isempty(regexpi(filename,'LOonly.raw','match'))
        continue
    end
    
    % Fetch number of measurement
    numberToken = regexpi(filename,'^([0123456789]*)-','tokens');
    number = str2double(cell2mat(numberToken{1}));
    dataStructLOonly(number).filename = filename;
    dataStructLOonly(number).number = number;
end
LOnumbers = cell2mat({dataStructLOonly.number});

for number = 1:size(dataStruct(:,2))
    filenameSIG = dataStruct(number).filename;
    %find adequate LOonly-file
    LOnumber = max(LOnumbers(LOnumbers<=number));    
    filenameLO = dataStructLOonly(LOnumber).filename;
    
    dispstat(['mainPrepareData number ' num2str(number)],'timestamp','keepthis',quiet);
    [ X, theta ] = mainPrepareData( filenameLO, filenameSIG );
    dispstat(['discretize data ' num2str(number)],'timestamp','keepthis',quiet);
    [ Xdis, thetadis ] = discretizeData( X, theta, 400 );
    dispstat(['compute and plot Expectations ' num2str(number)],'timestamp','keepthis',quiet);
    [ ~,~, ~, ~, expDelQ, expDelP, unc, ~, meanN ] =...
        computeExpectations2( Xdis, thetadis, strrep(num2str(filenameSIG),'-LOwithLO','') );
    dataStruct(number).meanN = meanN;
    dataStruct(number).unc = unc(1,1);
    dataStruct(number).delQ = expDelQ(1,1);
    dataStruct(number).delP = expDelP(1,1);   
end
meanNs = cell2mat({dataStruct.meanN});
uncs = cell2mat({dataStruct.unc});
delQs = cell2mat({dataStruct.delQ});
delPs = cell2mat({dataStruct.delP});

dispstat('plot uncertainties over mean photon number','timestamp','keepthis',quiet);
plot(meanNs,delQs,meanNs,delPs,meanNs,uncs, 'o');
legend( '$\Delta Q$','$\Delta P$', '$\Delta Q \cdot \Delta P$ ''Location','northeast');
axis ([0 max(meanNs)+1 0 2]);
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultTextInterpreter','latex');
xlabel('mean Photonnumber');
ylabel('Uncertainty [a. u.]');    
print('UncertaintiesOverPhotonnumber','-dpng');


end