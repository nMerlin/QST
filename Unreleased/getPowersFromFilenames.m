%% Variables
dataStruct = struct('filename',{},'I',{});

%% Create data overview
dispstat('','init');
dispstat('Checking filenames ...','timestamp','keepthis');
%Contents = dir('Quadratures');
Contents = dir('raw-data');
name = {Contents.name};

for iStruct =  1:length(Contents) 
    %get filename
    filename = cell2mat(name(iStruct));
    if  isempty(regexpi(filename,'LOwithDL','match'))...
        || not(isempty(regexpi(filename,'cfg','match')))...
        || not(isempty(regexpi(filename,'stamp','match')))
        continue
    end
    
%     if (isempty(regexpi(filename,'corrRemove-yes','match')))
%             continue
%     end
    
    
    dataStruct(iStruct).filename = filename;
    % get current or power
    %currentToken = regexpi(filename,'-([0123456789.]*)mA','tokens');
%     currentToken = regexpi(filename,'mW-([0123456789.]*)mW','tokens');
%      currentToken = regexpi(filename,'([0123456789,]*)mW-5mW','tokens');    
     %currentToken = regexpi(filename,'([0123456789,]*)mW','tokens'); 
    % currentToken = regexpi(filename,'mat([0123456789,]*)','tokens');
       
            %currentToken = regexpi(filename,'([0123456789,]*)mW','tokens');
    currentToken = regexpi(filename,'([0123456789,]*)mW-4mW','tokens');
     currentToken{1}=strrep(currentToken{1},',','.');
     dataStruct(iStruct).I = str2double(cell2mat(currentToken{1}));
     numberToken = regexpi(filename,'^([0123456789,]*)-','tokens'); 
     dataStruct(iStruct).number = str2double(cell2mat(numberToken{1}));
  
         
    
    
end % iStruct

Is = cell2mat({dataStruct.I});
numbers = cell2mat({dataStruct.number});
[~,I] = sort(numbers);
Is = Is(I);
save('powers.mat','numbers','Is');
