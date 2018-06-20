function T = series3Ch(varargin)
%SERIES3CH is used to batch process 3-Channel data already converted to mat
%
% Functionality:
%   Looks for *.mat files in the folder 'mat-data', computes multiple
%   quantities of interest, and saves them as a table. The
%   *.mat files should contain the vectors X1,X2 and X3 and a piezoSign
%   variable.
%
% Optional Input Arguments:
%   *** Postselection *** 'SavePostselection': Default is 'false'. Choose
%       'true' if you want to recompute and save the postselected variables
%       'O1', 'O2','O3','oTheta', 'selX', 'selTheta' and 'selParams' in a
%       separate file.
%   'SelectionParameters': Parameters for 'selectRegion'. Default is
%       'Type'='fullcircle' and 'Position'=[2.5,0.5]
%
%   *** Density Matrix ***
%   'RhoParams': Default is empty. Structure containing optional input
%       arguments for function 'computeDensityMatrix'.
%   'SaveRho': Default is 'false'. Choose true if you want the function to
%       compute the most likely density matrix in the selected data. If
%       such a matrix already exists, it's just loaded from memory.
%   'RewriteRho': Default is 'false'. Choose 'true' if you do not want to
%       read the density matrix from disk but want a fresh calculation.
%
%   *** Wigner Function ***
%   'SaveWigner': Default is 'false'. With 'true' the Wigner function will
%       be calculated. This also sets 'SaveRho' to 'true'.
%   'RewriteWigner': Default is 'false'. Choose 'true' if you do not want
%       to read the Wigner function from disk but a new calculation. This
%       also sets 'SaveRho' to 'true'.
%
%   *** Other ***
%   'SaveTheta': Default is 'false'. Choose true if you want the function
%       to update the *.mat file with the computed theta vector.
%   'NDisc': Default is 100. Number of bins to discretize theta for
%       computing expectation values with 'computeExpectations'.

%% Validate and parse input arguments
p = inputParser;
defaultFittedExpectations = false;
addParameter(p,'FittedExpectations',defaultFittedExpectations,@islogical);
defaultNDisc = 100;
addParameter(p,'NDisc',defaultNDisc,@isnumeric);
defaultRewriteRho = false;
addParameter(p,'RewriteRho',defaultRewriteRho,@islogical);
defaultRewriteWigner = false;
addParameter(p,'RewriteWigner',defaultRewriteWigner,@islogical);
defaultRhoParams = struct;
addParameter(p,'RhoParams',defaultRhoParams,@isstruct);
defaultSavePostselection = false;
addParameter(p,'SavePostselection',defaultSavePostselection,@islogical);
defaultSaveRho = false;
addParameter(p,'SaveRho',defaultSaveRho,@islogical);
defaultSaveTheta = false;
addParameter(p,'SaveTheta',defaultSaveTheta,@islogical);
defaultSaveWigner = false;
addParameter(p,'SaveWigner',defaultSaveWigner,@islogical);
defaultSelParams = struct('Type','fullcircle','Position',[2.5,0.5]);
addParameter(p,'SelectionParameters',defaultSelParams,@isstruct);
parse(p,varargin{:});
c = struct2cell(p.Results);
[fitexp,nDisc,rewriteRho,rewriteWigner,rhoParams,saveps,saverho, ...
    savetheta,saveWigner,selParams] = c{:};

% Dependencies among optional input arguments
if saveWigner || rewriteWigner
    saverho = true;
end

%% Discover *.mat files
filestruct = dir('mat-data/*.mat');
files = {filestruct.name};

%% Create folder 'post-data'
if ~exist('post-data','dir')
    mkdir('post-data')
end

%% Iterate through data files
quantities = struct; % Structure that will contain quantities of interest
tempsaveps = false; % Can change postselection saving behavior per case
selStr = selParamsToStr(selParams);
dispstat('','init','timestamp','keepthis',0);
for i = 1:length(files)
    %% Load data
    dispstat(['Processing ',files{i},' ...'],'timestamp','keepthis',0);
    clear X1 X2 X3 theta piezoSign
    clear O1 O2 O3 oTheta selX selTheta;
    clear rho WF;
    C = strsplit(files{i},'.');
    filename = C{1};
    if ~saveps
        try
            load(['post-data/',filename,'-',selStr,'.mat']);
        catch
            dispstat(['Could not find postselection file, ', ...
                'loading raw quadratures ...'],'timestamp','keepthis',0);
            load(['mat-data/',files{i}]);
            tempsaveps = true;
        end
    else
        load(['mat-data/',files{i}]);
    end

    %% Compute Phase and Postselected Variables
    if ~exist('selX','var') % run only if postselection file was not loaded
        if ~exist('theta','var') % run only if theta is not available
            try
                [theta,~] = computePhase(X1,X2,piezoSign);
            catch
                warning(['Problem using computePhase.', ...
                    ' Assigning random phase.']);
                theta = rand(size(X1,1)*size(X1,2),size(X1,3))*2*pi;
            end
        end
        [O1,O2,O3,oTheta] = selectOrthogonal(X2,X3,X1,theta,piezoSign);
        [selX,selTheta] = selectRegion(O1,O2,O3,oTheta,selParams);
        close all;
        
        % Compute photon numbers for each channel
        [n1,n2,n3] = nPhotons(X1,X2,X3);
        quantities.nX1(i) = n1;
        quantities.nX2(i) = n2;
        quantities.nX3(i) = n3;
        quantities.minVar(i) = compute3ChLimit(n1,n2,n3);
    end
    
    %% Get quantities of interest
    % Compute expectation values of postselected state by fitting
    if fitexp
        expectations = computeExpectationsFit(selX,selTheta);
        quantities.cohAmpl(i) = expectations.cohAmpl;
        quantities.fitMeanVar(i) = expectations.varX;
    end
    
    % Compute expectation values of postselected state
    [disSelX,disSelTheta]=discretizeTheta(selX,selTheta,nDisc);
    [expQ,expP,expQ2,expP2,delQ,delP,~,~,~,~] = ...
        computeExpectations2(disSelX,disSelTheta,'bla','Plot','hide');
    quantities.q1(i) = expQ(1);
    quantities.q21(i) = expQ2(1);
    quantities.p1(i) = expP(1);
    quantities.p21(i) = expP2(1);
    quantities.maxQ(i) = max(expQ);
    quantities.q2max(i) = max(expQ2);
    quantities.varQ(i) = delQ^2;
    quantities.varP(i) = delP^2;
    quantities.discAmpl(i) = expQ(round(nDisc/4));
    quantities.discMeanVar(i) = mean(expQ2(:)-(expQ(:)).^2);
    
    %% Save workspace variables (because recomputing them takes time)
    % Save Theta
    if savetheta
        save(['mat-data/',files{i}],'X1','X2','X3','theta','piezoSign');
    end
    
    % Save postselected variables
    if saveps || tempsaveps
        postFilename = ['post-data/',filename,'-',selStr,'.mat'];
        save(postFilename,'O1','O2','O3','oTheta', ...
            'selX','selTheta','selParams');
        tempsaveps = false;
        if exist('rho','var')
            save(postFilename,'rho','rhoParams','-append');
        end
        if exist('WF','var')
            save(postFilename,'WF','-append');
        end
    end
    
    %% Reconstruct the Quantum State
    % Compute the density matrix
    if (saverho && (~exist('rho','var'))) || rewriteRho
        rho = computeDensityMatrix(selX,selTheta,rhoParams);
        save(['post-data/',filename,'-',selStr,'.mat'], ...
            'rho','rhoParams','-append');
    end
    
    % Compute Wigner function
    if (saveWigner && ~exist('WF','var')) || rewriteWigner
        WF = mainWignerFromRho(rho);
        save(['post-data/',filename,'-',selStr,'.mat'], ...
            'WF','-append');
    end
end

%% Create and write table
% Load most recent table file 'yyyy-MM-dd-series3Ch.txt'
T = seriesRead3ChTable();
if isempty(T)
    try
        T = readtable('raw-data/Delays.txt');
    catch
        warning('There was no ''raw-data/Delays.txt''!');
    end
end

% Update table with new values by looping over 'quantities' variable
fields = fieldnames(quantities);
for iField = 1:numel(fields)
    T.(fields{iField}) = makecol(quantities.(fields{iField}));
end

% Remove '.mat' in filenames and write results to a new table file
[~,T.Filename]=cellfun(@fileparts,files','UniformOutput',false);
writetable(T,[datestr(date,'yyyy-mm-dd-'),'series3Ch-',selStr,'.txt']);

end
