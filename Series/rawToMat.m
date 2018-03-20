function rawToMat(varargin)
%RAWTOMAT Computing quadratures from *.raw files and saving them as *.mat
%
% This function should be run one folder above "raw-data" (housing *.raw)
% and will create a "mat-data" folder, if necessary. For each signal file
% there should be a 'LOonly' file with a lower number.
%
% Currently, only preparation of phase averaged 1-channel measurements
% supported.

%% Validate and parse input arguments
p = inputParser;
defaultPrepOpts = struct;
addParameter(p,'PrepOpts',defaultPrepOpts,@isstruct);
parse(p,varargin{:});
c = struct2cell(p.Results);
[prepopts] = c{:};

%% Create folder 'mat-data'
if ~exist('mat-data','dir')
    mkdir('mat-data')
end

%% Read file names and separate LO and Signal numbering
filestruct = dir('raw-data/*.raw');
filestring = strjoin({filestruct.name});
fLO = regexpi(filestring,'[^ ]*LOonly[^ ]*','match');
tokLO = regexpi(strjoin(fLO),'\<(\d*)-','tokens');
nLO = cellfun(@str2num,[tokLO{:}]);
fSig = regexpi(filestring,'[^ ]*-(?!LOonly)[^ ]*','match');
fSig = {fSig{~ismember(fSig,fLO)}}; % remove false positives
tokSig = regexpi(strjoin(fSig),'\<(\d*)-','tokens');
nSig = cellfun(@str2num,[tokSig{:}]);

%% Compute quadratures and save them
for i=1:length(fSig)
    % find index of corresponding LO (lower number)
    [~,iLO] = max(nLO(nLO<nSig(i)));
    filename = strsplit(fSig{i});
    dispstat(['Working on File ',filename{1}],'timestamp','keepthis',0);
    X = preparePhAvData(fLO{iLO},fSig{i},prepopts);
    parsave('mat-data', ...
        [datestr(date,'yyyy-mm-dd'),'-',filename{1},'.mat'],X);
end

% Leave a note about how the files were created
cd('mat-data');
fid = fopen('created-by-rawToMat.txt','a');
fclose(fid);

end

