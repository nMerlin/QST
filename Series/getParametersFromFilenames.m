function [filenames,numbers,Is]= getParametersFromFilenames(varargin)
% this function derives the parameter I and the number from filenames that are called like
% "number-I-4mW-LOwithDL-...". The files are found in a folder that can be
% specified. 

%% Validate and parse input arguments
p = inputParser;
defaultFolder = 'raw-data'; 
addParameter(p,'Folder',defaultFolder);
defaultParameter = 'power'; % which parameter was changed during the series
addParameter(p,'Parameter',defaultParameter);
defaultEnding = 'mat'; % ending of the filenames
addParameter(p,'Ending',defaultEnding);
parse(p,varargin{:});
c = struct2cell(p.Results);
[ending,folder,parameter] = c{:};

%% Create data overview
dispstat('','init');
dispstat('Checking filenames ...','timestamp','keepthis');
Contents = dir([folder '/*.' ending]);
name = {Contents.name};
filenames = name;

 dataStruct = struct('I',{});

for iStruct =  1:length(Contents) 
    %get filename
    filename = cell2mat(name(iStruct));
    if  isempty(regexpi(filename,'LOwithDL','match'))...
        || not(isempty(regexpi(filename,'cfg','match')))...
        || not(isempty(regexpi(filename,'stamp','match')))...
        || isempty(filename)
      continue
    end
    
%     if (isempty(regexpi(filename,'corrRemove-yes','match')))
%             continue
%     end
    
    dataStruct(iStruct).filename = filename;
    switch parameter
        case 'power'    
            %currentToken = regexpi(filename,'([0123456789,]*)mW','tokens');
            currentToken = regexpi(filename,'([0123456789,]*)mW-4mW','tokens');
             currentToken{1}=strrep(currentToken{1},',','.');
             dataStruct(iStruct).I = str2double(cell2mat(currentToken{1}));
            numberToken = regexpi(filename,'^([0123456789,]*)-','tokens'); 
             dataStruct(iStruct).number = str2double(cell2mat(numberToken{1}));
        case 'delay'
            delayToken = regexpi(filename,'([-0123456789,-]*)mm','tokens');
            delay = cell2mat(delayToken{1});
            numberToken = regexpi(delay,'^([0123456789,]*)-','tokens'); 
             dataStruct(iStruct).number = str2double(cell2mat(numberToken{1}));
            number = cell2mat(numberToken{1});
            delay = strrep(delay,[number '-'],'');
            delay = strrep(delay,',','.');
            delay = str2double(delay);
            c = 299792458; % in m/s
            delay = 2*delay/1000/c*10^12; %delay in ps   
            dataStruct(iStruct).I = delay;
        case 'no' 
            dataStruct(iStruct).I = 0;
            numberToken = regexpi(filename,'^([0123456789,]*)-','tokens'); 
             dataStruct(iStruct).number = str2double(cell2mat(numberToken{1}));             
        case 'position'
            delayToken = regexpi(filename,'([-0123456789,-]*)mm','tokens');
            delay = cell2mat(delayToken{1});
            numberToken = regexpi(delay,'^([0123456789,]*)-','tokens'); 
             dataStruct(iStruct).number = str2double(cell2mat(numberToken{1}));
            number = cell2mat(numberToken{1});
            delay = strrep(delay,[number '-'],'');
            delay = strrep(delay,',','.');
            delay = str2double(delay);
            dataStruct(iStruct).I = delay;
    end
     
end % iStruct

numbers = cell2mat({dataStruct.number});
Is = cell2mat({dataStruct.I});
[~,I] = sort(numbers);
Is = Is(I);
filenames = filenames(I);
save('parameters.mat','filenames','numbers','Is');
end
