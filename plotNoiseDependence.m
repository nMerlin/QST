function [ powerLO, deltaQ ] = plotNoiseDependence( varargin )
%DETSTANDARDTEST analyses the quadrature data of a detector standard test.
%   DETSTANDARDTEST('verbose'): Shows log output.
%   POWERLO: Processed LO powers.
%   DELTAQ: Calculated distribution widths.
%
%   See also LOAD8BITBINARY.

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

% Parameters & Variables
windowSize = 40; %Integrationwindow
fitThreshold = 5; %fitting a*sqrt(x) only for powers >= fitThreshold
outputFilename = 'shot-noise-plot.jpg'; %name of resulting saved plot
filetype = '-djpeg';

dataStruct = struct('filename',{},'powerLO',{},'deltaQ',{});

%%% Create data overview
dispstat('','init',quiet);
dispstat('Checking filenames ...','timestamp','keepthis',quiet);
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
dispstat('Calculating deltaQ ...','timestamp','keepthis',quiet);
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
set(0,'DefaultTextInterpreter','latex');
xlabel('P_{LO} [mW]');
ylabel('\Delta Q [a. u.]')
legend('Experimental $\Delta Q$ Data: $\sqrt{Var(Q)}$',strcat('Fit result ($P_{LO} \geq ',num2str(fitThreshold),'$ mW): $',num2str(sqrtFit.a),'*\sqrt{P_{LO}}$'),'Electronic Noise: $\Delta Q(0$ mW$)$','Location','northwest');
hold off;

% Saving figure
print(outputFilename,filetype);

dispstat('Finished!','timestamp',quiet);

end