function [ deltaQ, powerLO ] = plotShotNoise( varargin )
%PLOTSHOTNOISE analyses the quadrature data of a detector standard test.
%   PLOTSHOTNOISE('verbose'): Shows log output.
%   POWERLO: Processed LO powers.
%   DELTAQ: Calculated distribution widths.
%   The script assumes that the datafiles are in the 'raw-data' directory.
%   The filename-convention is '03-0.1mW-*.raw'.
%
%   See also LOAD8BITBINARY.
%
%   Optional Input Arguments:
%     'Verbose': Option 'quiet' mutes the command line output.

%% Validate and parse input arguments
p = inputParser;
defaultCompensate = false;
addParameter(p,'Compensate',defaultCompensate,@islogical);
defaultDutyCycle = 0;
addParameter(p,'DutyCycle',defaultDutyCycle,@isfloat);
defaultExclude = [];
addParameter(p,'Exclude',defaultExclude,@isvector);
defaultQuiet = 'notquiet';
addParameter(p,'Verbose',defaultQuiet,@isstr);
parse(p,varargin{:});
c = struct2cell(p.Results);
[compensate,dutycycle,exclude,quiet] = c{:};

if dutycycle>0
    intOpts.DutyCycle = dutycycle;
end

%% Constants
fitThreshold = 5; %fitting a*sqrt(x) only for powers >= fitThreshold
wavelength = 800e-9;
repetitionRate = 75.4e6;
planck = 6.626070040e-34;
lightSpeed = 299792458;
amperePerVolt = 1e-09;
powerConversion = wavelength/(planck*lightSpeed*repetitionRate)...
    /1000; %from mW to #LO photons per pulse
outputFilename = 'shot-noise-plot';
outputFiletype = '-djpeg';

dataStruct = struct('filename',{},'powerLO',{},'NLO',{},'deltaQ',{});

%% Create data overview
dispstat('','init',quiet);
dispstat('Checking filenames ...','timestamp','keepthis',quiet);
rawDataContents = dir('raw-data');
for name = {rawDataContents.name}
    % Loop only over *.raw files
    filename = cell2mat(name);
    if not(isempty(regexpi(filename,'.raw.','match'))) || ...
            isempty(regexpi(filename,'.raw','match'))
        continue
    end
    
    % Fetch number of measurement
    numberToken = regexpi(filename,'^([0123456789]*)-','tokens');
    number = str2double(cell2mat(numberToken{1}));
    dataStruct(number).filename = filename;
    if ismember(number,exclude)
        continue;
    end
    
    % Get LO power
    powerToken = regexpi(filename,'-([0123456789.]*)mW','tokens');
    dataStruct(number).powerLO = str2double(cell2mat(powerToken{1}));
    dataStruct(number).NLO = dataStruct(number).powerLO*powerConversion;
end
powerLO = cell2mat({dataStruct.powerLO});
[maxPowerLO,~] = max(powerLO);

% Calculate deltaQ=sqrt(var(Q)) of the quadrature values and get positions
% for integration from each LO-power
dispstat('Calculating deltaQ ...','timestamp','keepthis',quiet);
% Use locations of max power for all integrations
[data8bitMax,config,~]=load8BitBinary(dataStruct(end).filename,'dontsave');
[intOpts.Locations,~]=pointwiseVariance(data8bitMax,'MinPeakDistance',10);
X = computeQuadratures(data8bitMax,config,amperePerVolt,intOpts);

% Optional compensation of correlations
if compensate
    % Remove Offset
    X = bsxfun(@minus, X, mean(X));
    % Compensate correlations
    X = correlationCompensation(X);
end

dataStruct(end).deltaQ = sqrt(var(X(:)));

parfor number=1:(size(dataStruct,2)-1)
    if ~ismember(number,exclude)
    [data8bit,~,~]=load8BitBinary(dataStruct(number).filename,'dontsave');
    X = computeQuadratures(data8bit,config,amperePerVolt,intOpts);
    dataStruct(number).deltaQ = sqrt(var(X(:)));
    end
end

%% Create shot noise plot
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
loglog(powerLO,deltaQ,'o','Color','black','MarkerFaceColor','black');
hold on;
set(gca,'FontSize',14);
loglog(plotX,fitY);
loglog(plotX,electronicNoise,'Color','blue');
hold off;
axis tight;
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultTextInterpreter','latex');
xlabel('LO Power (mW)');
ylabel('\Delta Q (arb. units)')
legend('$\Delta Q = \sqrt{Var(Q_i)}$',...
    strcat('Fit Result ($P_{LO} \geq ',num2str(fitThreshold),...
    '$ mW): $',num2str(sqrtFit.a),'*\sqrt{P_{LO}}$'),...
    'Noise Background: $\Delta Q(0$ mW$)$','Location','northwest');
ax1 = gca; % current axes
ax1Pos = ax1.Position;
scale = 0.9;
ax1Pos(2) = ax1Pos(2)+(1-scale)/2*ax1Pos(4);
ax1Pos(4) = scale*ax1Pos(4);
set(ax1, 'Position', ax1Pos);
%ax1.Box = 'off';
% ax2 = axes('Position',ax1Pos,...
%     'Box','off',...
%     'XAxisLocation','top',...
%     'XScale','log',...
%     'XLim',[ax1.XLim(1)*powerConversion ax1.XLim(2)*powerConversion],...
%     'YAxisLocation','right',...
%     'YScale','log',...
%     'YLim',ax1.YLim,...
%     'Color','none','FontSize',14);
% ax2.XLabel.String = 'N_{LO}';

% Saving figure
%print([outputFilename,'.jpg'],outputFiletype);
%saveA5Landscape(outputFilename);
set(gcf,'Color','w');
ax1.XRuler.TickLength = 10;
ax1.YRuler.TickLength = 10;
formatFigA5;
export_fig shot-noise-plot.pdf

dispstat('Finished!','timestamp',quiet);

end