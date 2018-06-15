function [] = series3Ch(varargin)
%SERIES3CH is used to batch process 3-Channel data already converted to mat
%
% Functionality:
%   Looks for *.mat files in the folder 'mat-data', computes multiple
%   quantities of interest, and saves them as a table. The
%   *.mat files should contain the vectors X1,X2 and X3 and a piezoSign
%   variable.
%
% Optional Input Arguments:
%   'SaveTheta': Default is 'false'. Choose true if you want the function
%       to update the *.mat file with the computed theta vector.
%   'SavePostselection': Default is 'false'. Choose 'true' if you want the
%       function to save not only the raw quadratures, but also in a
%       separate file the postselected variables 'O1', 'O2', 'O3',
%       'oTheta', 'selX', 'selTheta' and 'selParams'.

%% Validate and parse input arguments
p = inputParser;
defaultSavePostselection = false;
addParameter(p,'SavePostselection',defaultSavePostselection,@islogical);
defaultSaveTheta = false;
addParameter(p,'SaveTheta',defaultSaveTheta,@islogical);
parse(p,varargin{:});
c = struct2cell(p.Results);
[saveps,savetheta] = c{:};

%% Discover *.mat files
filestruct = dir('mat-data/*.mat');
files = {filestruct.name};

%% Iterate through data files
quantities = struct; % Structure that will contain quantities of interest
tempsaveps = false; % Can change postselection saving behavior per case
dispstat('','init','timestamp','keepthis',0);
for i = 1:length(files)
    %% Load data
    dispstat(['Processing ',files{i},' ...'],'timestamp','keepthis',0);
    clear X1 X2 X3 theta piezoSign O1 O2 O3 oTheta selX selTheta selParams;
    if ~saveps
        try
            load(['mat-data/',name,'-postselection']);
        catch
            dispstat(['Could not find postselection file, ', ...
                'loading raw quadratures ...'],'timestamp','keepthis',0);
            load(['mat-data/',files{i}]);
            tempsaveps = true;
        end
    else
        load(['mat-data/',files{i}]);
    end
    [~,filename] = fileparts(files{i});

    %% Compute Phase and Postselect
    if ~exist('theta','var')
        try
            [theta,~] = computePhase(X1,X2,piezoSign);
        catch
            warning('Problem using computePhase. Assigning random phase.');
            theta = rand(size(X1,1)*size(X1,2),size(X1,3))*2*pi;
        end
    end
    [O1,O2,O3,oTheta] = selectOrthogonal(X2,X3,X1,theta,piezoSign);
    selParams.Type = 'fullcircle';
    selParams.Position = [2.5 0.5];
    [selX,selTheta,~,mVarX] = selectRegion(O1,O2,O3,oTheta,selParams, ...
        'Plot','show','Output','print','Filename',filename);
    close all;
    meanVarX(i) = mVarX;
    
    %% Get quantities of interest
    % Compute photon numbers for each channel
    [n1,n2,n3] = nPhotons(X1,X2,X3);
    quantities.nX1(i) = n1;
    quantities.nX2(i) = n2;
    quantities.nX3(i) = n3;
    
    % Compute expectation values of postselected state by fitting
    expectations = computeExpectationsFit(selX,selTheta);
    quantities.cohAmpl(i) = expectations.cohAmpl;
    
    % Compute expectation values of postselected state
    [selX,selTheta]=discretizeTheta(selX,selTheta,100);
    [expQ,expP,expQ2,expP2,delQ,delP,~,~,~,~] = ...
        computeExpectations2(selX,selTheta,'bla','Plot','hide');
    quantities.q1(i) = expQ(1);
    quantities.q21(i) = expQ2(1);
    quantities.p1(i) = expP(1);
    quantities.p21(i) = expP2(1);
    quantities.qmax(i) = max(expQ);
    quantities.q2max(i) = max(expQ2);
    quantities.varQ(i) = delQ^2;
    quantities.varP(i) = delP^2;
    
    %% Save high-cost workspace variables
    % Save Theta
    if savetheta
        save(['mat-data/',files{i}],'X1','X2','X3','theta','piezoSign');
    end
    
    % Save postselected variables
    if saveps || tempsaveps
        [~,name] = fileparts(files{i});
        save(['mat-data/',name,'-postselection'],'O1','O2','O3', ...
            'oTheta','selX','selTheta','selParams');
        tempsaveps = false;
    end
end

%% Create and write table
% Load most recent table file 'yyyy-MM-dd-series3Ch.txt' (easier way?)
filestruct = dir('*-series3CH.txt');
if ~isempty(filestruct)
    filestring = strjoin({filestruct.name});
    filedates = regexp(filestring,'([^ ]*)-series3Ch.txt','tokens');
    filedates = [filedates{:}]';
    filedates = datetime(filedates,'InputFormat','yyyy-MM-dd');
    filenames = {filestruct.name}';
    T = table(filedates,filenames);
    T = sortrows(T,'filedates');
    T = readtable(T.filenames{end});
end

% Update table with new values by looping over 'quantities' variable
fields = fieldnames(quantities);
for iField = 1:numel(fields)
    T.(fields{iField}) = makecol(quantities.(fields{iField}));
end

% Write results to a new table file
T.Filename = files';
T.minVar = compute3ChLimit(quantities.nX1,quantities.nX2,quantities.nX3)';
writetable(T,[datestr(date,'yyyy-mm-dd-'),'series3Ch']);

end