function plotDispersionsPositionSeries(varargin)
% this script makes spectrum plots for a series of dispersion measurements.
% The data should be in folder 'raw-data'.
% With 'Subtract', you can choose whether there are background measurements
% that should be subtracted. 
% ZeroPosition: The position of the left edge of sample in mm. 

%% Validate and parse input arguments
parser = inputParser;
defaultSubtract = 'no'; %
addParameter(parser,'Subtract',defaultSubtract);
defaultPlottype = 'lin'; 
addParameter(parser,'Plottype',defaultPlottype);
defaultZeroPosition = 3.67; 
addParameter(parser,'ZeroPosition',defaultZeroPosition,@isnumeric);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[plottype,subtract,zeroPos] = c{:};

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
    E0 = plotDispersion(filenameSIG, filenameBG,'Subtract',subtract,'Plottype',plottype);
    dataStruct(number).E0 = E0;
       
end

position = cell2mat({dataStruct.Position});
position = position - zeroPos;
E0 = cell2mat({dataStruct.E0});

%% write them in excel table
T = struct2table(dataStruct);
writetable(T,'positionSeries.xls');

%% make plot
plot(position,E0);
xlabel('sample position relative to left edge (mm)');
ylabel('minimum energy E_{LP}(k_{||}=0)(eV)');
graphicsSettings;
grid();
savefig('positionSeries.fig');
print('positionSeries.png','-dpng','-r300');

end