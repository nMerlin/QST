function T = seriesRead3ChTable()
%SERIESREADTABLE3CH Load most recent table file 'yyyy-MM-dd-series3Ch.txt'
% in current folder.

filestruct = dir('*-series3Ch.txt');
if ~isempty(filestruct)
    filestring = strjoin({filestruct.name});
    filedates = regexp(filestring,'([^ ]*)-series3Ch.txt','tokens');
    filedates = [filedates{:}]';
    filedates = datetime(filedates,'InputFormat','yyyy-MM-dd');
    filenames = {filestruct.name}';
    T = table(filedates,filenames);
    T = sortrows(T,'filedates');
    T = readtable(T.filenames{end});
else
    warning('There is no table *-series3Ch.txt available!');
    T = [];
end

end

