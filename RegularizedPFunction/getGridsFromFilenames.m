function [xGrid,phiGrid] = getGridsFromFilenames(directory)
%gets xGrid and phiGrid from the filenames of archived pattern functions. 

Contents = dir([directory '/*.mat']);
name = {Contents.name};
[xPreliminary, phiPreliminary] = deal(zeros(length(Contents),1));
for i = 1:length(Contents)
    filename = cell2mat(name(i));
    xToken = regexpi(filename,'x-([-0123456789.]*)-','tokens');
    xPreliminary(i) = str2double(cell2mat(xToken{1}));
    phiToken = regexpi(filename,'phi-([-0123456789.]*).mat','tokens');
    phiPreliminary(i) = str2double(cell2mat(phiToken{1}));
end
xGrid= unique(xPreliminary)';
phiGrid = unique(phiPreliminary)';

end


