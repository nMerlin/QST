function T = seriesRead3ChTable(selParams)
%SERIESREADTABLE3CH Load most recent table file 'yyyy-MM-dd-series3Ch.txt'
% in current folder. Output is empty if no table was detected.
%
% Usage:
%   T = seriesRead3ChTable(selParams);
%
% Input Arguments:
%   selParams: Structure with selection parameters for 'selectRegion'.

selStr = '';
if nargin > 0
    selStr = selParamsToStr(selParams);
end
filestruct = dir(['*-series3Ch-',selStr,'.txt']);
if ~isempty(filestruct)
    filestring = strjoin({filestruct.name});
    filedates = regexp(filestring,['([^ ]*)','-series3Ch-',selStr] ...
        ,'tokens');
    filedates = [filedates{:}]';
    filedates = datetime(filedates,'InputFormat','yyyy-MM-dd');
    filenames = {filestruct.name}';
    T = table(filedates,filenames);
    T = sortrows(T,'filedates');
    T = readtable(T.filenames{end});
else
    warning(['There is no table *-series3Ch-',selStr,'.txt available!']);
    T = [];
end

end

