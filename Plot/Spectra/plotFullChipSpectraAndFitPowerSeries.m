function plotFullChipSpectraAndFitPowerSeries(varargin)
% this script makes spectrum plots for a series of spectrum measurements.
% The data should be in folder 'raw-data'.
% With 'Subtract', you can choose whether there are background measurements
% that should be subtracted. 
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
defaultSubtract = 'yes'; %
addParameter(parser,'Subtract',defaultSubtract);
defaultXLim = []; %
addParameter(parser,'XLim',defaultXLim,@isnumeric);
defaultXUnit = 'nm';
addParameter(parser,'XUnit',defaultXUnit);
defaultXrange = [];
addParameter(parser,'Xrange',defaultXrange);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[fitoption,intp,save,subtract,xLim,xRange,xUnit] = c{:};

%% Create data overview
dataStruct = struct('filename',{},'number',{}, 'Power',{}, 'Max', {}, 'Integrated',{});
dataStructBackground = struct('filename',{},'number',{});

rawDataContents = dir('raw-data');

%% Signal-files
for name = {rawDataContents.name}
    % Loop only over Signal files
    filename = cell2mat(name);
%     if not(isempty(regexpi(filename,'background_1.txt','match')))...
%             || isempty(regexpi(filename,'.txt','match'))
%         continue
%     end
    
    if not(isempty(regexpi(filename,'background','match')))...
            || isempty(regexpi(filename,'.csv','match'))
        continue
    end
    
    % Fetch number of measurement
    numberToken = regexpi(filename,'^([0123456789]*)-','tokens');
    number = str2double(cell2mat(numberToken{1}));
    dataStruct(number).filename = filename;
    
    % Fetch excitation power
    powerToken = regexpi(filename,'-([0123456789.]*)mW','tokens');
    %powerToken = regexpi(filename,'OPO-([0123456789.]*)mW','tokens');
    %powerToken = regexpi(filename,'MIRA-([0123456789.]*)mW','tokens');
    power = str2double(cell2mat(powerToken{1}));
    dataStruct(number).Power = power;
    
end

%% Background-files
for name = {rawDataContents.name}
    % Loop only over background files
    filename = cell2mat(name);
    if isempty(regexpi(filename,'background.txt','match'))
        continue
    end
    
    % Fetch number of measurement
    numberToken = regexpi(filename,'^([0123456789]*)-','tokens');
    number = str2double(cell2mat(numberToken{1}));
    dataStructBackground(number).filename = filename;
    dataStructBackground(number).number = number;
end
BGnumbers = cell2mat({dataStructBackground.number});

%% process the data
for number = 1:size(dataStruct,2)

    filenameSIG = dataStruct(number).filename;
    if isempty(filenameSIG)
        continue
    end    
   
    %find adequate BG-file, if this option is chosen
    if strcmp(subtract, 'yes')
        BGnumber = min(BGnumbers(BGnumbers>=number)); %background was measured after signal
        filenameBG = dataStructBackground(BGnumber).filename;
        [Max,integratedInt] = plotFullChipSpectrumAndFit( filenameSIG, filenameBG,...
            'Subtract',subtract,'Interpolate',intp,'Fit',fitoption,'Save',save,'XLim',xLim,'XRange',xRange,'XUnit',xUnit); 
    else
        [Max,integratedInt] = plotFullChipSpectrumAndFit( filenameSIG, filenameSIG,...
            'Subtract',subtract,'Interpolate',intp,'Fit',fitoption,'Save',save,'XLim',xLim,'XRange',xRange,'XUnit',xUnit); 
    end
    dataStruct(number).Max = Max;
    dataStruct(number).Integrated = integratedInt;
       
end

powers = cell2mat({dataStruct.Power});
maxs = cell2mat({dataStruct.Max});
Integrateds = cell2mat({dataStruct.Integrated});

%% write them in excel table
T = struct2table(dataStruct);
writetable(T,'powerSeries.xls');

%% make plots
loglog(powers,maxs,'o','LineWidth',2,'DisplayName','maximum peak Intensity');
l = legend('Location','northwest');
l.FontSize = 22;
title('PowerSeries-max');
ylabel('Counts');
xlabel('Excitation Power (mW)');
graphicsSettings;
grid();
savefig('PowerSeries-max.fig');
print('PowerSeries-max.png','-dpng','-r300');

loglog(powers,Integrateds,'o','LineWidth',2,'DisplayName','Int. Intensity');
l = legend('Location','northwest');
l.FontSize = 22;
title('PowerSeries-Integrated');
ylabel('Counts');
xlabel('Excitation Power (mW)');
graphicsSettings;
grid();
savefig('PowerSeries-Integrated.fig');
print('PowerSeries-Integrated.png','-dpng','-r300');
clf();

end