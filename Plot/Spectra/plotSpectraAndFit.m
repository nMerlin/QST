function plotSpectraAndFit(varargin)
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
defaultXLim = 0.5; %
addParameter(parser,'XLim',defaultXLim,@isnumeric);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[fitoption,intp,save,subtract,xLim] = c{:};

%% Create data overview
dataStruct = struct('filename',{},'number',{}, 'Max', {}, ...
    'peak', {}, 'FWHM', {}, 'Q', {});
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
        [Max,integratedInt, peak, FWHM, Q] = plotSpectrumAndFit( filenameSIG, filenameBG,...
            'Subtract',subtract,'Interpolate',intp,'Fit',fitoption,'Save',save,'XLim',xLim); 
    else
        [Max,integratedInt, peak, FWHM, Q] = plotSpectrumAndFit( filenameSIG, filenameSIG,...
            'Subtract',subtract,'Interpolate',intp,'Fit',fitoption,'Save',save,'XLim',xLim); 
    end
       
end



end