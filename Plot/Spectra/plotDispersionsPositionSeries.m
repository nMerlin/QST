function plotDispersionsPositionSeries(varargin)
% this script makes spectrum plots for a series of dispersion measurements.
% The data should be in folder 'raw-data'.
% With 'Subtract', you can choose whether there are background measurements
% that should be subtracted. 
% ZeroPosition: The position of the left edge of sample in mm. 

%% Validate and parse input arguments
parser = inputParser;


%% Validate and parse input arguments
parser = inputParser;
defaultSubtract = 'no'; %
addParameter(parser,'Subtract',defaultSubtract);
defaultPlottype = 'lin'; 
addParameter(parser,'Plottype',defaultPlottype);
defaultPlotMode = 'no'; 
addParameter(parser,'PlotMode',defaultPlotMode); % if yes, the selected mode is indicated in the plot. 
defaultMinY = 200; 
addParameter(parser,'minY',defaultMinY,@isnumeric);
defaultXAperture = 330'; 
addParameter(parser,'xAperture',defaultXAperture,@isnumeric);
defaultModeK = [-0.66 0.66]; 
addParameter(parser,'ModeK',defaultModeK,@isnumeric); % k range in which the mode is selected
defaultModeE = [1.611 1.6113]; 
addParameter(parser,'ModeE',defaultModeE,@isnumeric); % Energy range in which the mode is selected
defaultZoomE = [1.605 1.625]; 
addParameter(parser,'ZoomE',defaultZoomE,@isnumeric); % Energy range in which the dispersion is plotted
defaultFit = 'yes'; %
addParameter(parser,'Fit',defaultFit);
defaultE0Old = 1.6108; 
addParameter(parser,'E0Old',defaultE0Old,@isnumeric);
defaultaOld = 4.4973e-7; 
addParameter(parser,'aOld',defaultaOld,@isnumeric);
defaultY0Old = 217; 
addParameter(parser,'y0Old',defaultY0Old,@isnumeric);
defaultGetTime = 'no'; %
addParameter(parser,'GetTime',defaultGetTime); % get exposure time from filename
defaultLabel = '-'; %
addParameter(parser,'Label',defaultLabel); % use only files with specific label e.g. polarisation
defaultZeroPosition = 3.67; 
addParameter(parser,'ZeroPosition',defaultZeroPosition,@isnumeric);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[aOld,E0Old,fitoption,getTime,label,minY,modeE,modeK,plotMode,plottype,subtract,xAperture,y0Old,zeroPos,zoomE] = c{:};

%% Sample char.
R = 9.5; %Rabi splitting in meV
EX = 1.6195; % Exciton energy in eV

%% Create data overview
dataStruct = struct('filename',{},'number',{}, 'Position',{}, 'E0', {});
dataStructBackground = struct('filename',{},'number',{});

rawDataContents = dir('raw-data');

%% Signal-files
for name = {rawDataContents.name}
    % Loop only over Signal files
    filename = cell2mat(name);
    if not(isempty(regexpi(filename,'background_1.txt','match')))...
            || isempty(regexpi(filename,'.txt','match'))
        continue
    end
    
    % Fetch number of measurement
    numberToken = regexpi(filename,'^([0123456789]*)-','tokens');
    number = str2double(cell2mat(numberToken{1}));
    dataStruct(number).filename = filename;
    
    % Fetch sample position
    positionToken = regexpi(filename,'-([0123456789.]*)mm','tokens');
    position = str2double(cell2mat(positionToken{1}));
    dataStruct(number).Position = position;
    
end

%% Background-files
for name = {rawDataContents.name}
    % Loop only over background files
    filename = cell2mat(name);
    if isempty(regexpi(filename,'background_1.txt','match'))
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
    else
        filenameBG = filenameSIG;
    end
    E0 = plotDispersion(filenameSIG, filenameBG,'Subtract',subtract,...
        'Plottype',plottype,'minY',minY,'xAperture',xAperture,'ModeE',modeE,...
        'ModeK',modeK,'ZoomE',zoomE,'Fit',fitoption,'aOld',aOld,'E0Old',E0Old,'y0Old',y0Old,'PlotMode',plotMode);
    dataStruct(number).E0 = E0;
       
end

position = cell2mat({dataStruct.Position});
position = position - zeroPos;
E0 = cell2mat({dataStruct.E0});
[detuningMeV, detuning2g0] = computeDetuning(E0,R,EX); 

%% write them in excel table
T = struct2table(dataStruct);
writetable(T,'positionSeries.xls');

%% make plot of Energy vs position
plot(position,E0);
xlabel('sample position relative to left edge (mm)');
ylabel('minimum energy E_{LP}(k_{||}=0)(eV)');
graphicsSettings;
grid();
savefig('positionSeries-EnergyVsPosition.fig');
print('positionSeries-EnergyVsPosition.png','-dpng','-r300');

%% make plot of Energy vs detuning (meV)
plot(detuningMeV,E0);
xlabel('detuning (meV)');
ylabel('minimum energy E_{LP}(k_{||}=0)(eV)');
graphicsSettings;
grid();
savefig('positionSeries-EnergyVsDetuningMeV.fig');
print('positionSeries-EnergyVsDetuningMeV.png','-dpng','-r300');

%% make plot of Energy vs detuning (2g0)
plot(detuning2g0,E0);
xlabel('detuning (2 g_{0})');
ylabel('minimum energy E_{LP}(k_{||}=0)(eV)');
graphicsSettings;
grid();
savefig('positionSeries-EnergyVsDetuning2g0.fig');
print('positionSeries-EnergyVsDetuning2g0.png','-dpng','-r300');

%% make plot of detuning (meV) vs position
plot(position,detuningMeV);
xlabel('sample position relative to left edge (mm)');
ylabel('detuning (meV)');
graphicsSettings;
grid();
savefig('positionSeries-detuningMeVvsPosition.fig');
print('positionSeries-detuningMeVvsPosition.png','-dpng','-r300');

%% make plot of detuning (2g0) vs position
plot(position,detuningMeV);
xlabel('sample position relative to left edge (mm)');
ylabel('detuning (2 g_{0})');
graphicsSettings;
grid();
savefig('positionSeries-detuning2g0vsPosition.fig');
print('positionSeries-detuning2g0vsPosition.png','-dpng','-r300');

end