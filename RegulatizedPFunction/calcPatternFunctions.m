function [] = calcPatternFunctions(directory)
%for each set of X parameters, this calculates the pattern functions and
%saves them. It gets R from folder name. 

%% Normierung???

%get R from folder name
RToken = regexpi(directory,'-R-([0123456789.]*)','tokens');
R = cell2mat(RToken{1});
R = str2double(R);

%get filenames 
Contents = dir([directory '/*XParameter.mat']);
name = {Contents.name};

for iStruct =  1:length(Contents) 
    %get filename
    filename = [directory '\' cell2mat(name(iStruct))];
    load(filename,'XGrid');
    % compute the pattern function from the X parameters, see Jans notes, equ. 14 and 15
    pattern = patternFunction(XGrid,R);
    newFilename = strrep(filename,'XParameter','Pattern');
    save(newFilename, 'pattern');
end
   