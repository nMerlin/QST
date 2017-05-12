function [ deltaQ, powerLO ] = plotShotNoise( varargin )
%PLOTSHOTNOISE analyses the quadrature data of a detector standard test.
%   PLOTSHOTNOISE('verbose'): Shows log output.
%   POWERLO: Processed LO powers.
%   DELTAQ: Calculated distribution widths.
%   The script assumes that the datafiles are in the 'raw-data' directory.
%   The filename-convention is '03-0.1mW-*.raw'.
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
wavelength = 800e-9;
repetitionRate = 75.4e6;
planck = 6.626070040e-34;
lightSpeed = 299792458;
powerConversion = wavelength/(planck*lightSpeed*repetitionRate)...
    /1000; %from mW to #LO photons per pulse
outputFilename = 'shot-noise-plot.jpg';
outputFiletype = '-djpeg';

dataStruct = struct('filename',{},'powerLO',{},'NLO',{},'deltaQ',{});

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
    dataStruct(number).NLO = dataStruct(number).powerLO*powerConversion;
end


powerLO = cell2mat({dataStruct.powerLO});
[maxPowerLO,~] = max(powerLO);

% Calculate deltaQ=sqrt(var(Q)) of the quadrature values;  Get positions for integration from each LO-power
dispstat('Calculating deltaQ ...','timestamp','keepthis',quiet);
%For 0mW, the locations of 1 mW are used.
[data8bit0,~,~] = load8BitBinary(dataStruct(1).filename,'dontsave');  
[data8bit1,~,~] = load8BitBinary(dataStruct(2).filename,'dontsave');
[locs,~] = pointwiseVariance(data8bit1);
[~,X]=correlation(0,data8bit0,locs);
dataStruct(1).deltaQ = sqrt(var(X(:)));

parfor number=2:size(dataStruct,2)
    [data8bit,~,~] = load8BitBinary(dataStruct(number).filename,'dontsave');  
    [locs,~] = pointwiseVariance(data8bit);
    [~,X]=correlation(0,data8bit,locs);
    dataStruct(number).deltaQ = sqrt(var(X(:)));
end

%%% Create shot noise plot
% Process data
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

dispstat('Plotting ...','timestamp','keepthis',quiet);
% Plotting
close all;
figure;
loglog(powerLO,deltaQ,'o');
hold on;
loglog(plotX,fitY);
loglog(plotX,electronicNoise);
hold off;
axis tight;
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultTextInterpreter','latex');
xlabel('P_{LO} [mW]');
ylabel('\Delta Q [a. u.]')
legend('Experimental $\Delta Q$ Data: $\sqrt{Var(Q)}$',...
    strcat('Fit Result ($P_{LO} \geq ',num2str(fitThreshold),...
    '$ mW): $',num2str(sqrtFit.a),'*\sqrt{P_{LO}}$'),...
    'Electronic Noise: $\Delta Q(0$ mW$)$','Location','northwest');
ax1 = gca; % current axes
ax1Pos = ax1.Position;
scale = 0.9;
ax1Pos(2) = ax1Pos(2)+(1-scale)/2*ax1Pos(4);
ax1Pos(4) = scale*ax1Pos(4);
set(ax1, 'Position', ax1Pos);
ax1.Box = 'off';
ax2 = axes('Position',ax1Pos,...
    'Box','off',...
    'XAxisLocation','top',...
    'XScale','log',...
    'XLim',[ax1.XLim(1)*powerConversion ax1.XLim(2)*powerConversion],...
    'YAxisLocation','right',...
    'YScale','log',...
    'YLim',ax1.YLim,...
    'Color','none');
ax2.XLabel.String = 'N_{LO}';

% Saving figure
print(outputFilename,outputFiletype);

dispstat('Finished!','timestamp',quiet);

end