function plotSpectraAndFit
% to do: background subtraction optional
% this script makes spectrum plots for a series of spectrum measurements.
% The data should be in folder 'raw-data'.

%% Create data overview
dataStruct = struct('filename',{});
dataStructBackground = struct('filename',{},'number',{}, 'Max', {}, ...
    'peak', {}, 'FWHM', {}, 'Q', {});

rawDataContents = dir('raw-data');

%% Signal-files
for name = {rawDataContents.name}
    % Loop only over Signal files
    filename = cell2mat(name);
    if not(isempty(regexpi(filename,'background.dat','match')))...
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
   
    %find adequate BG-file
%     BGnumber = min(BGnumbers(BGnumbers>=number)); %background was measured after signal
%     filenameBG = dataStructBackground(BGnumber).filename;

    [Max, peak, FWHM, Q] = plotSpectrumAndFit( filenameSIG);
       
end



end