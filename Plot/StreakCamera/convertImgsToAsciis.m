function [] =  convertImgsToAsciis
% for a series of streak measurements located in folder 'raw-data'.
% Converts .img files to .dat files (ascii).


%% Create data overview
rawDataContents = dir('raw-data');

cd('raw-data');

for name = {rawDataContents.name}
    % Loop only over Signal files
    filename = cell2mat(name);
    if isempty(regexpi(filename,'.img','match'))
        continue
    end
    
    loadImg(filename); 
    
end

cd('..');


end