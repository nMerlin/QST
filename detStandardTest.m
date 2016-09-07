function [ dataStruct ] = detStandardTest( )
%DETSTANDARDTEST analyses the shot noise behavior from a detector standard
%test.
%   Details of data format etc.

% Parameters
windowSize = 40; %Integrationwindow
fitThreshold = 5; %fitting a*sqrt(x) only for powers >=fitThreshold

% Variables
dataStruct = struct('filename',{},'powerLO',{},'deltaQ',{});

%%% Create an data overview
dispstat('','init');
dispstat('Checking filenames ...','timestamp','keepthis');
rawDataContents = dir('raw-data');
for name = {rawDataContents.name}
    % Loop only over *.raw files
    filename = cell2mat(name);
    if not(isempty(regexpi(filename,'.raw.','match'))) || isempty(regexpi(filename,'.raw','match'))
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

% Get positions for integration from highest LO-power
powerLO = cell2mat({dataStruct.powerLO});
[maxPowerLO,index] = max(powerLO);
[data8bit,~,~] = load8BitBinary(dataStruct(index).filename,'dontsave');
[locs,~] = pointwiseVariance(data8bit);

% Calculate deltaQ=sqrt(var(Q)) of the quadrature values
dispstat('Calculating deltaQ ...','timestamp','keepthis');
parfor number=1:size(dataStruct,2)
    [data8bit,~,~] = load8BitBinary(dataStruct(number).filename,'dontsave');
    [~,X]=correlation(0,data8bit,locs,windowSize);
    dataStruct(number).deltaQ = sqrt(var(X(:)));
end

%%% Create shot noise plot
plotX = 0.1:0.1:maxPowerLO;
deltaQ = cell2mat({dataStruct.deltaQ});
[~,index] = min(powerLO);
electronicNoise = ones(size(plotX))*dataStruct(index).deltaQ;

% Fitting
fitIndices = find(powerLO>=fitThreshold);
fitPowerLO = transpose(powerLO(fitIndices));
fitDeltaQ = transpose(deltaQ(fitIndices));
ft = fittype('a*sqrt(x)');
sqrtFit = fit(fitPowerLO,fitDeltaQ,ft,'StartPoint',5);
fitY = sqrtFit(plotX);

% Plotting
loglog(powerLO,deltaQ,'o');
hold on;
loglog(plotX,fitY);
loglog(plotX,electronicNoise);
set(0,'DefaultLegendInterpreter','latex');
legend('1','2','$\Delta Q(0$ mW$)$','Location','northwest');
hold off;

dispstat('Finished!','timestamp');

end