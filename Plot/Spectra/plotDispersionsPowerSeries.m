function plotDispersionsPowerSeries(varargin)
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
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[plottype,subtract] = c{:};

%% Create data overview
dataStruct = struct('filename',{},'number',{}, 'Power',{}, 'E0', {}, 'Emax', {}, 'Intmax', {}, 'SumInt', {});
dataStructBackground = struct('filename',{},'number',{});

rawDataContents = dir('raw-data');

%% Signal-files
for name = {rawDataContents.name}
    % Loop only over Signal files
    filename = cell2mat(name);
    if not(isempty(regexpi(filename,'-Background.csv','match')))...
            || isempty(regexpi(filename,'.csv','match'))
        continue
    end
    
    % Fetch number of measurement
    numberToken = regexpi(filename,'^([0123456789]*)-','tokens');
    number = str2double(cell2mat(numberToken{1}));
    dataStruct(number).filename = filename;
    
    % Fetch excitation power
    powerToken = regexpi(filename,'-([0123456789.]*)mW','tokens');
    power = str2double(cell2mat(powerToken{1}));
    dataStruct(number).Power = power;
    
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
    [E0,Emax,Intmax,SumInt] = plotDispersion(filenameSIG, filenameBG,'Subtract',subtract,'Plottype',plottype);
    dataStruct(number).E0 = E0;
    dataStruct(number).Emax = Emax; 
    dataStruct(number).Intmax = Intmax; 
    dataStruct(number).SumInt = SumInt; 
end

power = cell2mat({dataStruct.Power});
E0 = cell2mat({dataStruct.E0});
Emax = cell2mat({dataStruct.Emax});
Intmax = cell2mat({dataStruct.Intmax});
SumInt = cell2mat({dataStruct.SumInt});

%% write them in excel table
T = struct2table(dataStruct);
writetable(T,'powerSeries.xls');

%% make plots
plot(power,E0);
xlabel('Excitation Power (mW)');
ylabel('minimum energy E_{LP}(k_{||}=0)(eV)');
graphicsSettings;
grid();
savefig('powerSeries-E0.fig');
print('powerSeries-E0.png','-dpng','-r300');
clf();

plot(power,Emax,'o');
xlabel('Excitation Power (mW)');
ylabel('Energy of maximum (eV)');
graphicsSettings;
grid();
savefig('powerSeries-Emax.fig');
print('powerSeries-Emax.png','-dpng','-r300');
clf();

%% constants
h = 6.62607015e-34;
c0 = 299792458;
e0 = 1.602176634e-19;
lambda = h*c0/e0*1./Emax *10^9;
plot(power,lambda,'o');
xlabel('Excitation Power (mW)');
ylabel('wavelength of maximum (nm)');
graphicsSettings;
grid();
savefig('powerSeries-wavelength.fig');
print('powerSeries-wavelength.png','-dpng','-r300');
clf();

loglog(power,Intmax,'o');
xlabel('Excitation Power (mW)');
ylabel('maximum integrated Intensity (a.u.)');
graphicsSettings;
grid();
savefig('powerSeries-Intmax.fig');
print('powerSeries-Intmax.png','-dpng','-r300');
clf();

end