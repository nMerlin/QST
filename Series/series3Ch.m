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
%   'ChannelAssignment': default is [1,2,3]; This sets which channel is target, 
%       which is the postselection channel that is piezo modulated fast for 
%       the phase computation and which is the last postselection channel, 
%       according to [target,ps_piezo_fast,ps_piezo_slow]
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
%       to update the *.mat file with the computed theta vector, and with 
%       the orhtogonal variables O1,O2,O3,oTheta.
%   'RecomputeTheta': Default is 'false'. Choose 'true' if you do not want
%       to read the theta from disk but a new calculation, which is not saved. 
%   'SaveOrth': Default is 'false'. Choose true if you want the function
%       to update the *.mat file with the orthogonal variables O1,O2,O3,oTheta.
%   'RecomputeOrth': Default is 'false'. Choose 'true' if you do not want
%       to read the orthogonal variables O1,O2,O3,oTheta from disk but a new 
%       calculation, which is not saved. 
%   'NDisc': Default is 100. Number of bins to discretize theta for
%       computing expectation values with 'computeExpectations'.
%   'Parameter': Default is 'delay'. Then, delay is retrieved
%       from the file name with format xx-yymm, where yy is the delay.
%       choose 'power' when power was changed during the series and get
%       this frome filename.
%   'RemoveModulation': Default is 'false'. Choose true if only data should
%       be selected where the photon number lies in a certain range 
%   'range': the range for selecting only photon numbers within this range
%       (it takes nPsFast into regard)
%   'VaryAPS': Default is 'false'. choose the postselection radius according to 
%       photon numbers, while keeping the conditional radius as fixed parameter
%   'CorrRemove': Only uses those mat-file that have 'corrRemove-yes' in
%       the filename
%   'Period': number of periods expected in one piezo segment. Important for phase computation. 
%
%   *** Developer Only ***
%   'FileRange': Default is [] and loops over all found files. In any other
%       case, the loop is only executed over the given file indices. This
%       is intended for manual parallelization only.

%% Validate and parse input arguments
p = inputParser;
defaultFileRange = []; % be careful
addParameter(p,'FileRange',defaultFileRange,@isvector);
defaultFileNumber = 0; % be careful
addParameter(p,'FileNumber',defaultFileNumber,@isnumeric);
defaultFittedExpectations = false;
addParameter(p,'FittedExpectations',defaultFittedExpectations,@islogical);
defaultNDisc = 100;
addParameter(p,'NDisc',defaultNDisc,@isnumeric);
defaultPeriod = 2; % number of periods expected in one piezo segment. Important for phase computation. 
addParameter(p,'Period',defaultPeriod,@isnumeric);
defaultRewriteRho = false;
addParameter(p,'RewriteRho',defaultRewriteRho,@islogical);
defaultRewriteWigner = false;
addParameter(p,'RewriteWigner',defaultRewriteWigner,@islogical);
defaultRecomputeTheta = false;
addParameter(p,'RecomputeTheta',defaultRecomputeTheta,@islogical);
defaultRecomputeOrth = false;
addParameter(p,'RecomputeOrth',defaultRecomputeOrth,@islogical);
defaultRhoParams = struct;
addParameter(p,'RhoParams',defaultRhoParams,@isstruct);
defaultSavePostselection = false;
addParameter(p,'SavePostselection',defaultSavePostselection,@islogical);
defaultSaveRho = false;
addParameter(p,'SaveRho',defaultSaveRho,@islogical);
defaultSaveTheta = false;
addParameter(p,'SaveTheta',defaultSaveTheta,@islogical);
defaultSaveOrth = false;
addParameter(p,'SaveOrth',defaultSaveOrth,@islogical);
defaultSaveWigner = false;
addParameter(p,'SaveWigner',defaultSaveWigner,@islogical);
defaultParameter = 'delay';
addParameter(p,'Parameter',defaultParameter,@isstr);
defaultRemoveModulation = false;
addParameter(p,'RemoveModulation',defaultRemoveModulation,@islogical);
defaultVaryAPS = false;
addParameter(p,'VaryAPS',defaultVaryAPS,@islogical);
defaultCorrRemove = 'yes';
addParameter(p,'CorrRemove',defaultCorrRemove,@isstr);
defaultSelParams = struct('Type','fullcircle','Position',[2.5,0.5]);
addParameter(p,'SelectionParameters',defaultSelParams,@isstruct);
defaultRange = [0 10];
addParameter(p,'Range',defaultRange,@isvector);
defaultChannelAssignment = [3,1,2]; %[target,ps_piezo_fast,ps_piezo_slow]
addParameter(p,'ChannelAssignment',defaultChannelAssignment,@isvector);
parse(p,varargin{:});
c = struct2cell(p.Results);
[chAssign,corrRemove,filenumber,filerange,fitexp,nDisc,parameter,periodsPerSeg,range,recomputeOrth,recomputeTheta,remMod,rewriteRho,rewriteWigner,rhoParams, ...
    saveOrth,saveps,saverho,savetheta,saveWigner,selParams,varyAPS] = c{:};

% Dependencies among optional input arguments
if saveWigner || rewriteWigner
    saverho = true;
end

%% Discover *.mat files
if strcmp(corrRemove,'yes')
        filestruct = dir('mat-data/*CorrRemove-yes.mat');
    else
        filestruct = dir('mat-data/*.mat'); %% will discover all files !
end


files = {filestruct.name};

%% Create folder 'post-data'
if ~exist([pwd 'post-data'],'dir')
    mkdir('post-data')
end

%% Iterate through data files
quantities = struct; % Structure that will contain quantities of interest
tempsaveps = false; % Can change postselection saving behavior per case
selStr = selParamsToStr(selParams);
dispstat('','init','timestamp','keepthis',0);
if isempty(filerange)
    filerange = 1:length(files);
end
for ii = 1:length(filerange)
    i = filerange(ii);
    %% Load data
    C = strsplit(files{i},'.');
    filename = C{1};
    dispstat(['Loading ',files{i},' ...'],'timestamp',0);
    clear X1 X2 X3 theta piezoSign
    clear O1 O2 O3 oTheta selX selTheta;
    clear rho WF;
    postFilename =  ['post-data/',filename,'-',selStr,'-remMod-',...
            num2str(remMod),'-range-',num2str(range),'-varyAPS-',num2str(varyAPS),'.mat'];
    
    switch parameter
        case 'power'    
            powerToken =  regexpi(filename,'([0123456789,]*)mW-4mW','tokens');
            power = cell2mat(powerToken{1});
            power = strrep(power,',','.');
            power = str2double(power);
            quantities.Power(i) = power;      
        case 'delay'
        % get delay from the file name with format xx-yymm, where xx is the
        % filenumber and yy is the delay which can also be negative and start
        % with a minus sign        
            delayToken = regexpi(filename,'([-0123456789,-]*)mm','tokens');
            delay = cell2mat(delayToken{1});
            numberToken = regexpi(delay,'^([0123456789,]*)-','tokens'); 
            number = cell2mat(numberToken{1});
            delay = strrep(delay,[number '-'],'');
            delay = strrep(delay,',','.');
            delay = str2double(delay);
            quantities.DelayMm(i) = delay;  % in mm
            c = 299792458; % in m/s
            delay = 2*delay/1000/c*10^12; %delay in ps
            quantities.Delay(i) = delay;      
    end
    
    if ~saveps && ~recomputeTheta
        try
            load(postFilename);
        catch
            dispstat(['Could not find ',postFilename, ...
                ' loading raw quadratures ...'],'timestamp','keepthis',0);
            load(['mat-data/',files{i}]);
            tempsaveps = true;
        end
    else
        load(['mat-data/',files{i}]);
    end
    
    if exist('X1','var')
        quadratures = zeros([size(X1) 3]);
        quadratures(:,:,:,1) = X1;
        quadratures(:,:,:,2) = X2;
        quadratures(:,:,:,3) = X3;
        Xtg = quadratures(:,:,:,chAssign(1));  % this sets X1 to be the target channels etc
        XpsFast = quadratures(:,:,:,chAssign(2));
        XpsSlow = quadratures(:,:,:,chAssign(3));
        clear('quadratures'); 
    end 
    
    %% Compute Phase and Postselected Variables
    if ~exist('selX','var') % run only if postselection file was not loaded
        if (~exist('theta','var')) || recomputeTheta % run only if theta is not available or should be rewritten
            try
                [theta,~] = computePhase(Xtg,XpsFast,piezoSign,'Period',periodsPerSeg);
            catch
                warning(['Problem using computePhase.', ...
                    ' Assigning random phase.']);
                theta = rand(size(Xtg,1)*size(Xtg,2),size(Xtg,3))*2*pi;
            end
            
        end
              
        if (~exist('O1','var')) || recomputeOrth 
            [O1,O2,O3,oTheta,iOrth] = selectOrthogonal(XpsFast,XpsSlow,Xtg,theta,piezoSign);
        end
        
        if remMod 
                nTgVec = photonNumberVector(Xtg);
                nPsFastVec = photonNumberVector(XpsFast);
                nPsSlowVec = photonNumberVector(XpsSlow);
                %nTgVec = nTgVec(iOrth);
                nPsFastVec = nPsFastVec(iOrth);
                %nPsSlowVec = nPsSlowVec(iOrth);
                iSel = find(nPsFastVec >= min(range) & nPsFastVec <= max(range));
                O1 = O1(iSel);
                O2 = O2(iSel);
                O3 = O3(iSel);
                oTheta = oTheta(iSel);                
               % [O3,O1,O2,oTheta] = removeNBelowLimit(O3,O1,O2,oTheta,nTgVec,nPsFastVec,nPsSlowVec,range);
%             end
        end
        
        % Compute photon numbers for each channel
        [nTg,nPsFast,nPsSlow] = nPhotons(Xtg,XpsFast,XpsSlow);       
        
        selParamsUse = selParams;
        if varyAPS    %choose the postselection radius according to photon numbers
           selParamsUse.Position(1) = selParams.Position(1) * (1+nPsFast+nPsSlow)/sqrt(2*nTg*(nPsFast+nPsSlow));
           quantities.APS(i) = selParamsUse.Position(1);
        end
               
        [selX,selTheta] = selectRegion(O1,O2,O3,oTheta,selParamsUse);%,'Plot','show','Filename',[filename '-assessTheta']
        close all;
        
    end
    
    %% Get quantities of interest
    % Compute expectation values of postselected state by fitting
    if (exist('nTg','var'))
        quantities.nTg(i) = nTg;
        quantities.nPsFast(i) = nPsFast;
        quantities.nPsSlow(i) = nPsSlow;
        quantities.minVar(i) = compute3ChLimit(nPsFast,nPsSlow,nTg);
    end
    
    if fitexp
        expectations = computeExpectationsFit(selX,selTheta);
        quantities.cohAmpl(i) = expectations.cohAmpl;
        quantities.fitMeanVar(i) = expectations.varX;
    end
    
    % Compute simple expectation values of postselected state
    [disSelX,disSelTheta]=discretizeTheta(selX,selTheta,nDisc);
    [expQ,expP,expQ2,expP2,~,~,~,~,~,~] = ...
        computeExpectations2(disSelX,disSelTheta,'bla','Plot','hide');
    quantities.q1(i) = expQ(1);
    quantities.q21(i) = expQ2(1);
    quantities.p1(i) = expP(1);
    quantities.p21(i) = expP2(1);
    quantities.maxQ(i) = max(expQ);
    quantities.q2max(i) = max(expQ2);
    quantities.discAmpl(i) = mean(expQ(round(nDisc/4)-2:round(nDisc/4)+2));
    quantities.discMeanVar(i) = mean(expQ2(:)-(expQ(:)).^2);
    varVector =  expQ2(:)-(expQ(:)).^2; 
    quantities.varQ(i) = mean(mean([varVector(round(nDisc/4)-2:round(nDisc/4)+5) varVector(round(nDisc*3/4)-2:round(nDisc*3/4)+5)]));
    quantities.varP(i) = mean(varVector(round(nDisc/2)-2:round(nDisc/2)+5));
     %quantities.discN(i) = 0.5 * (mean(expQ2)+mean(expP2) - 1);
    
    % g2 estimation
    [g2values,nValues] = deal(zeros(10,1));
    for iG2=1:10
        try
            uniformX = seriesUniformSampling(selX,selTheta,'NBins',100);
        catch
            warning(['Problem with uniformSampling.', ...
                'Use X without uniform sampling for g2.']);
            uniformX = selX;
        end    
        try
            [g2values(iG2),nValues(iG2)] = g2(uniformX,length(uniformX));
        catch
            warning('not enough X to compute g2. Set g2 = 2.');
            g2values(iG2) = 2;
            uniformX=uniformX-mean(uniformX);
            nValues(iG2) = mean(uniformX.^2)-0.5;
        end
    end
    quantities.g2(i) = mean(g2values);
    quantities.g2std(i) = std(g2values);
    quantities.n(i) = mean(nValues);
    meang2 = mean(g2values);
    meanVar = mean(expQ2(:)-(expQ(:)).^2);
    quantities.discN(i) = mean(nValues);
    
    %% Save workspace variables (because recomputing them takes time)
    % Save Theta
    if savetheta 
            save(['mat-data/',files{i}],'theta','-append');
    end
    
    if saveOrth
            save(['mat-data/',files{i}],'O1','O2','O3','oTheta','iOrth','-append');
    end
    
    % Save postselected variables
    if saveps || tempsaveps
        save(postFilename, ...
            'selX','selTheta','selParams','meang2','meanVar','nTg','nPsFast','nPsSlow');      
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
        save(postFilename,'rho','rhoParams','-append');
    end
    
    % Compute Wigner function
    if (saveWigner && ~exist('WF','var')) || rewriteWigner
        WF = mainWignerFromRho(rho);
        save(postFilename,'WF','-append');
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
        T = cell2table(files','VariableNames',{'Filename'});
    end
end

% Update table with new values by looping over 'quantities' variable
fields = fieldnames(quantities);
for iField = 1:numel(fields)
    T.(fields{iField}) = makecol(quantities.(fields{iField}));
end

% Remove '.mat' in filenames and write results to a new table file
files = strrep(files,',','.');
[~,T.Filename]=cellfun(@fileparts,files','UniformOutput',false);
writetable(T,[datestr(date,'yyyy-mm-dd-'),'series3Ch-',selStr,'-remMod-',...
    num2str(remMod),'-range-',num2str(range),'-varyAPS-',num2str(varyAPS),'.txt']);

end
