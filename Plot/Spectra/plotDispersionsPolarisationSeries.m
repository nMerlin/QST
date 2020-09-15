function plotDispersionsPolarisationSeries(varargin)
% this script makes spectrum plots for a series of dispersion measurements.
% The data should be in folder 'raw-data'.
% With 'Subtract', you can choose whether there are background measurements
% that should be subtracted. 

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
defaultFit = 'yes'; %'useOld' to use old parameters
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
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[aOld,E0Old,fitoption,getTime,label,minY,modeE,modeK,plotMode,plottype,subtract,xAperture,y0Old,zoomE] = c{:};

%% Sample char. for M3396
R = 9.5; %Rabi splitting in meV
EX = 1.6195; % Exciton energy in eV
%% constants
h = 6.62607015e-34;
c0 = 299792458;
e0 = 1.602176634e-19;

%% Create data overview
dataStruct = struct('filename',{},'number',{}, 'Polarisation',{}, 'time', {}, 'E0', {}, 'lambda',{}, 'Emax', {}, 'modeInt', {}, 'SumInt', {});
dataStructBackground = struct('filename',{},'number',{});

rawDataContents = dir('raw-data');

%% Signal-files
for name = {rawDataContents.name}
    % Loop only over Signal files
    filename = cell2mat(name);
    if not(isempty(regexpi(filename,'-background.csv','match')))...
             || not(isempty(regexpi(filename,'raw','match')))...
            || isempty(regexpi(filename,'.csv','match'))...
            || isempty(regexpi(filename,label,'match'))
        continue
    end
    
    % Fetch number of measurement
    numberToken = regexpi(filename,'^([0123456789]*)-','tokens');
    number = str2double(cell2mat(numberToken{1}));
    dataStruct(number).filename = filename;
    
    % Fetch polarisation
    polToken = regexpi(filename,'l2-([0123456789.]*)','tokens');
    polarisation = str2double(cell2mat(polToken{1}));
    dataStruct(number).Polarisation = polarisation*2; %position of halfwaveplate * 2 is polarisation 
    
    if strcmp(getTime,'yes')
        timeToken = regexpi(filename,'-([0123456789.]*)ms','tokens');
        time = str2double(cell2mat(timeToken{1}));
        dataStruct(number).time = time;
    else
        dataStruct(number).time = 1;
    end
    
end

%% Background-files
for name = {rawDataContents.name}
    % Loop only over background files
    filename = cell2mat(name);
    if isempty(regexpi(filename,'background','match'))
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
        BGnumber = min(BGnumbers(BGnumbers<=number)); %background was measured before signal
        filenameBG = dataStructBackground(BGnumber).filename;
    else
        filenameBG = filenameSIG;
    end
    [E0,~,~,Emax,modeInt,SumInt,FWHM,FWHMerror,peakPosition] = plotDispersion(filenameSIG, filenameBG,'Subtract',subtract,...
        'Plottype',plottype,'minY',minY,'xAperture',xAperture,'ModeE',modeE,...
        'ModeK',modeK,'ZoomE',zoomE,'Fit',fitoption,'aOld',aOld,'E0Old',E0Old,'y0Old',y0Old,'PlotMode',plotMode,'R',R,'EX',EX);
    lambda = h*c0/e0*1./E0 *10^9;
    dataStruct(number).E0 = E0;
    dataStruct(number).lambda = lambda;
    dataStruct(number).Emax = Emax; 
    dataStruct(number).FWHM =FWHM; 
    dataStruct(number).FWHMerror = FWHMerror; 
    dataStruct(number).peakPosition = peakPosition; 
    if strcmp(getTime,'yes')
        dataStruct(number).modeInt = modeInt/dataStruct(number).time;
        dataStruct(number).SumInt = SumInt/dataStruct(number).time;
    else
        dataStruct(number).modeInt = modeInt; 
        dataStruct(number).SumInt = SumInt;
    end
    
end

polarisation = cell2mat({dataStruct.Polarisation});
E0 = cell2mat({dataStruct.E0});
lambda = cell2mat({dataStruct.lambda});
Emax = cell2mat({dataStruct.Emax});
modeInt = cell2mat({dataStruct.modeInt});
SumInt = cell2mat({dataStruct.SumInt});
FWHM = cell2mat({dataStruct.FWHM});
FWHMerror = cell2mat({dataStruct.FWHMerror});
peakPosition = cell2mat({dataStruct.peakPosition});

%% write them in excel table
T = struct2table(dataStruct);
writetable(T,['powerSeries' label '.xls']);

%% make plots
[polarisation,Index] = sort(polarisation);
SumInt = SumInt(Index);
polar(polarisation*pi/180,SumInt,'-o');
xlabel('linear polarisation angle');
if strcmp(getTime,'yes')
   ylabel('overall integrated Intensity (counts/ms)'); 
else
    ylabel('overall integrated Intensity (a.u.)');
end
graphicsSettings;
savefig(['polarisation-polarPlot-SumInt' label '.fig']);
print(['polarisation-polarPlot-SumInt' label '.png'],'-dpng','-r300');
clf();

plot(polarisation,SumInt,'-o');
xlabel('linear polarisation angle');
if strcmp(getTime,'yes')
   ylabel('overall integrated Intensity (counts/ms)'); 
else
    ylabel('overall integrated Intensity (a.u.)');
end
graphicsSettings;
savefig(['polarisation-linearPlot-SumInt' label '.fig']);
print(['polarisation-linearPlot-SumInt' label '.png'],'-dpng','-r300');
clf();

modeInt = modeInt(Index);
polar(polarisation*pi/180,modeInt,'-o');
xlabel('linear polarisation angle');
if strcmp(getTime,'yes')
   ylabel('mode integrated Intensity (counts/ms)'); 
else
    ylabel('mode integrated Intensity (a.u.)');
end
graphicsSettings;
savefig(['polarisation-polarPlot-modeInt' label '.fig']);
print(['polarisation-polarPlot-modeInt' label '.png'],'-dpng','-r300');
clf();

plot(polarisation,modeInt,'-o');
xlabel('linear polarisation angle');
if strcmp(getTime,'yes')
   ylabel('mode integrated Intensity (counts/ms)'); 
else
    ylabel('mode integrated Intensity (a.u.)');
end
graphicsSettings;
savefig(['polarisation-linearPlot-modeInt' label '.fig']);
print(['polarisation-linearPlot-modeInt' label '.png'],'-dpng','-r300');
clf();
end