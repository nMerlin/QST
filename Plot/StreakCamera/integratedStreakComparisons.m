function integratedStreakComparisons(varargin)
% for a series of integrated streak measurements located in folder 'integrated-data' as .dat
% files. It plots the streak data of converted an unconverted light over each other and
% also the difference using the function integratedStreakComparisons. 

%% Validate and parse input arguments
parser = inputParser;
defaultSubtract = 'yes'; 
addParameter(parser,'Subtract',defaultSubtract);
defaultPlottype = 'lin'; 
addParameter(parser,'Plottype',defaultPlottype);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[plottype,subtract] = c{:};

%% Create data overview
dataStructConv = struct('filename',{},'Power',{},'decaytime',{},'decaytimeError',{}, 'Max', {}, 'Sum', {});
dataStructUnc = struct('filename',{},'Power',{},'decaytime',{},'decaytimeError',{}, 'Max', {}, 'Sum', {});

intDataContents = dir('integrated-data');

%% Converted files
number = 1;
for name = {intDataContents.name}
    % Loop only over Signal files
    filename = cell2mat(name);
    if isempty(regexpi(filename,'-converted-','match'))...
            || isempty(regexpi(filename,'.mat','match'))
        continue
    end
    
    dataStructConv(number).filename = filename;
    
     % Fetch excitation power
    %powerToken = regexpi(filename,'-([0123456789.]*)mW','tokens'); %mira
    powerToken = regexpi(filename,'MIRA-([0123456789.]*)mW','tokens'); %mira
    power = str2double(cell2mat(powerToken{1}));
    dataStructConv(number).Power = power; 
    
    number = number + 1;
end

%% Unconverted files
number = 1;
for name = {intDataContents.name}
    
    % Loop only over Signal files
    filename = cell2mat(name);
    if isempty(regexpi(filename,'-unconverted-','match'))...
            || isempty(regexpi(filename,'.mat','match'))
        continue
    end
    
    dataStructUnc(number).filename = filename;
    
     % Fetch excitation power
    %powerToken = regexpi(filename,'-([0123456789.]*)mW','tokens'); %mira
    powerToken = regexpi(filename,'MIRA-([0123456789.]*)mW','tokens'); %mira
    power = str2double(cell2mat(powerToken{1}));
    dataStructUnc(number).Power = power; 
    
    number = number + 1;
end

%power array
powersUnc = cell2mat({dataStructUnc.Power});


%% process the data
for number = 1:size(dataStructConv,2)
    
    filenameConv = dataStructConv(number).filename;
    if isempty(filenameConv)
        continue
    end    
   
    %find Unc file with the same power
    powerConv = dataStructConv(number).Power;
    numberUnc = find(powersUnc == powerConv);
    filenameUnc = dataStructUnc(numberUnc).filename;
    if isempty(filenameUnc)
        continue
    end
    
    %% make the comparison
    integratedStreakComparison(filenameConv,filenameUnc);
     
%     dataStruct(number).decaytime = decaytime; 
%     dataStruct(number).decayError = decaytimeError;
%     dataStruct(number).Max = Max;
%     dataStruct(number).Sum = Sum;
end


end


























