function plotSpectraAndFitDepletionPowerSeries(varargin)
% this script makes spectrum plots for a series of spectrum measurements.
% The data should be in folder 'raw-data'.
%       'Fit': Decide with 'yes', if the data should be
%       fitted with a gaussian.
%       'Interpolate': Decide with 'yes', if the data should be
%       interpolated with a spline. This interpolation will also be used
%       for the fit.
%       'XLim': limits for x-axis around the peak.

%% Validate and parse input arguments
parser = inputParser;
defaultFit = 'yes'; %
addParameter(parser,'Fit',defaultFit);
defaultSave = 'yes'; %
addParameter(parser,'Save',defaultSave);
defaultInterpolate = 'yes';
addParameter(parser,'Interpolate',defaultInterpolate);
defaultXLim = 5; %
addParameter(parser,'XLim',defaultXLim,@isnumeric);
defaultXUnit = 'nm';
addParameter(parser,'XUnit',defaultXUnit);
defaultXrange = [];
addParameter(parser,'Xrange',defaultXrange);
defaultToken = '-'; %a part of the filename, so only these files are processed
addParameter(parser,'Token',defaultToken);
defaultPowerAfterFactor = 0; %if this is zero, the OPO power after the waveguide is taken from filename.
% Otherwise, it is computed from the power before the waveguide. 
addParameter(parser,'PowerAfterFactor',defaultPowerAfterFactor);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[fitoption,intp,powerAfterFactor,saveOpt,token,xLim,xRange,xUnit] = c{:};

%% Create data overview
dataStruct = struct('filename',{},'number',{}, 'Power',{}, 'Max', {}, 'Integrated',{},...
    'peak', {}, 'FWHM', {}, 'Q', {});
rawDataContents = dir('raw-data');

%% Signal-files
for name = {rawDataContents.name}
    % Loop only over Signal files
    filename = cell2mat(name);
    
    if not(isempty(regexpi(filename,'background','match')))...
            || isempty(regexpi(filename,'.ssm','match'))...
            || isempty(regexpi(filename,token,'match'))
                    %|| isempty(regexpi(filename,'.csv','match'))...
        continue
    end
    
    % Fetch number of measurement
    numberToken = regexpi(filename,'^([0123456789]*)-','tokens');
    number = str2double(cell2mat(numberToken{1}));
    dataStruct(number).filename = filename;
    
    % Fetch OPO power before and after waveguide 
    powerToken1 = regexpi(filename,'OPO-([0123456789.]*)mW','tokens');
    powerBefore = str2double(cell2mat(powerToken1{1}));
    dataStruct(number).PowerBefore = powerBefore;
    if powerAfterFactor == 0
        powerToken2 = regexpi(filename,'mW-([0123456789.]*)mW','tokens');
        powerAfter = str2double(cell2mat(powerToken2{1}));
    else 
        powerAfter = powerAfterFactor*powerBefore;
    end
    dataStruct(number).PowerAfter = powerAfter;
    
    % process spectrum 
     [Max,integratedInt, peak, FWHM, Q] = plotSpectrumAndFit( filename, filename,...
        'Subtract','no','Interpolate',intp,'Fit',fitoption,'Save',saveOpt,'XLim',xLim,'XRange',xRange,'XUnit',xUnit);
    
    % decide whether OPO open or not
    if not(isempty(regexpi(filename,'-o-','match')))
        dataStruct(number).MaxO = Max;
        dataStruct(number).IntegratedO = integratedInt;
        dataStruct(number).peakO = peak;
    elseif not(isempty(regexpi(filename,'-c-','match')))
        dataStruct(number).MaxC = Max;
        dataStruct(number).IntegratedC = integratedInt;
        dataStruct(number).peakC = peak;
    end
      
end %name

powersBefore = cell2mat({dataStruct.PowerBefore});
powersAfter = cell2mat({dataStruct.PowerAfter});
maxsO = cell2mat({dataStruct.MaxO});
IntegratedsO = cell2mat({dataStruct.IntegratedO});
peaksO = cell2mat({dataStruct.peakO});
maxsC = cell2mat({dataStruct.MaxC});
IntegratedsC = cell2mat({dataStruct.IntegratedC});
peaksC = cell2mat({dataStruct.peakC});
depletionMax = maxsO./maxsC;
depletionInt = IntegratedsO./IntegratedsC;
efficiencyMax = 1 -  depletionMax;
efficiencyInt = 1 -  depletionInt;

%% write them in excel table
T = struct2table(dataStruct);
writetable(T,['OPOSeries-' token '.xls']);
%save matlab file 
save(['OPOseries-data-' token '.mat'],'powersBefore','powersAfter','maxsO','IntegratedsO','peaksO',...
    'maxsC','IntegratedsC','peaksC','depletionMax','depletionInt','efficiencyMax','efficiencyInt');

%% make plots Depletion max
plot(powersBefore,depletionMax*100,'o','markerFaceColor','b');
ax1=gca;
graphicsSettings;
set(ax1,'xlim',[0 max(powersBefore)],'Box','on');
ax2=axes('Position',ax1.Position,'XAxisLocation','top','YAxisLocation',...
    'right','color','none','YColor','none','XColor','r','linewidth',2,...
        'FontSize',20,'FontName','Arial');
set(ax2,'xlim',[0 max(powersAfter)]);
xlabel(ax1,'Pump power before waveguide (mW)');
xlabel(ax2,'Pump power after waveguide (mW)');
ylabel(ax1,'Depletion (%)');
savefig(['OPOSeries' token '-depletionMax.fig']);
print(['OPOSeries' token '-depletionMax.png'],'-dpng','-r300');
hold off; clf();

%% make plots Depletion Int
plot(powersBefore,depletionInt*100,'o');
l = legend('from integrated intensity','Location','northwest');
ylabel('Depletion (%)');
xlabel('Pump power (mW)');
graphicsSettings;
savefig(['OPOSeries' token '-depletionInt.fig']);
print(['OPOSeries' token '-depletionInt.png'],'-dpng','-r300');

%% make plots Efficiency Max
plot(powersBefore,efficiencyMax*100,'o','markerFaceColor','b');
ax1=gca;
graphicsSettings;
set(ax1,'xlim',[0 max(powersBefore)],'Box','on');
ax2=axes('Position',ax1.Position,'XAxisLocation','top','YAxisLocation',...
    'right','color','none','YColor','none','XColor','r','linewidth',2,...
        'FontSize',20,'FontName','Arial');
set(ax2,'xlim',[0 max(powersAfter)]);
xlabel(ax1,'Pump power before waveguide (mW)');
xlabel(ax2,'Pump power after waveguide (mW)');
ylabel(ax1,'Efficiency (%)');
savefig(['OPOSeries' token '-efficiencyMax.fig']);
print(['OPOSeries' token '-efficiencyMax.png'],'-dpng','-r300');
hold off; clf();

%% make plots Efficiency Int
plot(powersBefore,efficiencyInt*100,'o');
l = legend('from integrated intensity.png','Location','northwest');
ylabel('Efficiency (%)');
xlabel('Pump power (mW)');
graphicsSettings;
savefig(['OPOSeries' token '-efficiencyInt.fig']);
print(['OPOSeries' token '-efficiencyInt.png'],'-dpng','-r300');

%% plot maxima
plot(powersBefore,maxsO,'o',powersBefore,maxsC,'o');
l = legend('max with pump','max without pump', 'Location','northwest');
ylabel('Counts');
xlabel('Pump power (mW)');
graphicsSettings;
savefig(['OPOSeries' token '-max.fig']);
print(['OPOSeries' token '-max.png'],'-dpng','-r300');

%% plot integrated Intensity
plot(powersBefore,IntegratedsO,'o',powersBefore,IntegratedsC,'o');
l = legend('integrated with pump','integrated without pump','Location','northwest');
ylabel('Counts');
xlabel('Pump power (mW)');
graphicsSettings;
savefig(['OPOSeries' token '-Integrated.fig']);
print(['OPOSeries' token '-Integrated.png'],'-dpng','-r300');
clf();

%% plot peak wavelength
plot(powersBefore,peaksO,'o',powersBefore,peaksC,'o');
l = legend('peak position with pump','peak position without pump','Location','northwest');
ylabel('Wavelength (nm)');
xlabel('Pump power (mW)');
graphicsSettings;
savefig(['OPOSeries' token '-peaks.fig']);
print(['OPOSeries' token '-peaks.png'],'-dpng','-r300');


end