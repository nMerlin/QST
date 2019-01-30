function plotSpectraAndFitPolDegree(numberVector)
% to do: background subtraction optional
% to do: background subtraction optional
% this script makes spectrum plots for a series of spectrum measurements.
% For different polarisations. 
% The data should be in folder 'raw-data'.
% NUMBERVECTOR: the numbers of the datafiles that should be used.

%% Create data overview
dataStruct = struct('filename',{}, 'lambda2', {}, 'modeInt', {}, ...
    'modePeak', {}, 'FWHM', {}, 'Q', {}, 'statusMax',{});
dataStructBackground = struct('filename',{},'number',{});

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
    if number < numberVector(1) || number > numberVector(end)
        continue
    end
    
    %fetch position of lambda/2
    lambda2Token = regexpi(filename,'-([0123456789]*)degree','tokens');
    lambda2 = str2double(cell2mat(lambda2Token{1}));
    dataStruct(number).lambda2 = lambda2;
end

%% Background-files
% for name = {rawDataContents.name}
%     % Loop only over background files
%     filename = cell2mat(name);
%     if isempty(regexpi(filename,'background.txt','match'))
%         continue
%     end
%     
%     % Fetch number of measurement
%     numberToken = regexpi(filename,'^([0123456789]*)-','tokens');
%     number = str2double(cell2mat(numberToken{1}));
%     dataStructBackground(number).filename = filename;
%     dataStructBackground(number).number = number;
% end
% BGnumbers = cell2mat({dataStructBackground.number});

%% process the data
modePeak = 0;

for number = numberVector

    filenameSIG = dataStruct(number).filename;
    if isempty(filenameSIG)
        continue
    end    
   
    %find adequate BG-file
%     BGnumber = min(BGnumbers(BGnumbers>=number)); %background was measured after signal
%     filenameBG = dataStructBackground(BGnumber).filename;

%This assumes we start with a maximum measurement.
    [modeInt, peak, FWHM, Q] = plotSpectrumAndFit( filenameSIG);
    if not(isempty(regexpi(filenameSIG,'max','match'))) 
        modePeak = peak;
        dataStruct(number).status = 1; 
    end
    
      % for minimum position of lambda/2, we need the intensitiy of 
        %the main mode from the maxmimum position of lambda/2
    if not(isempty(regexpi(filenameSIG,'min','match'))) 
        dataStruct(number).status = 0; 
        cd('raw-data');
        data = textread(filenameSIG);
        w = data(:,1); % wavelength
        Int = data(:,2); % Intensity
        %subtract background
        bg = min(Int);
        Int = Int-bg;
        modeInt = Int(w==modePeak);
        cd('..');      
    end
    dataStruct(number).modeInt = modeInt;  
    dataStruct(number).modePeak = modePeak; 
    dataStruct(number).FWHM = FWHM;
    dataStruct(number).Q = Q; 
end

lambda2s = cell2mat({dataStruct.lambda2});
modeInts = cell2mat({dataStruct.modeInt});
modePeaks = cell2mat({dataStruct.modePeak});
FWHMs = cell2mat({dataStruct.FWHM});
Qs = cell2mat({dataStruct.Q});
statusMaxs = cell2mat({dataStruct.status});

maxima = modeInts(statusMaxs==1);
minima = modeInts(statusMaxs==0);
Qs = Qs(statusMaxs==1);
meanMax = mean(maxima);
meanMin = mean(minima);
stdMax = std(maxima);
stdMin = std(minima);
polDeg = (meanMax-meanMin)/(meanMax+meanMin);

%write table
filename = ['PolarisationDegree' num2str(numberVector(1)) '-' num2str(numberVector(end)) '.mat'];
save(filename, 'lambda2s', 'modeInts', 'modePeaks','FWHMs','Qs','statusMaxs','maxima','minima','polDeg','meanMax','meanMin');

end