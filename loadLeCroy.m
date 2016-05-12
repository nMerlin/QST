function [ data ] = loadLeCroy( regex )
%LOADLECROY reads datafiles that match the regular expression REGEX and
% loads them into the DATA matrix.

filenames = dir(regex);
data = [];

for i=1:length(flucfiles)
    data = [data, csvread(flucfiles(i).name,5)];
end

end

