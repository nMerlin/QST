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
nLO = cellfun(@str2num,regexpi(strjoin(fLO),'\<\d*','match'));
fSig = regexpi(filestring,'[^ ]*-(?!LOonly)[^ ]*','match');
nSig = cellfun(@str2num,regexpi(strjoin(fSig),'\<\d*','match'));

%% Compute quadratures and save them
for i=1:length(fSig)
    % find index of corresponding LO (lower number)
    iSig = nSig(i);
    [~,iLO] = max(nLO(nLO<iSig));
    X = preparePhAvData(fLO{iLO},fSig{iSig},prepopts);
    cd('mat-data');
    filename = strsplit(fSig{iSig});
    save([datestr(date,'yyyy-mm-dd'),'-',filename{1},'.mat'],'X');
    cd('..');
end

end

