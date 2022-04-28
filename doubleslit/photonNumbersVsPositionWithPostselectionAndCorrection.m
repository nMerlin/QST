function [] = photonNumbersVsPositionWithPostselectionAndCorrection(varargin)

%% Validate and parse input arguments
p = inputParser;
defaultFolder = 'mat-data'; % which folder contains the data
addParameter(p,'Folder',defaultFolder);
defaultChannelAssignment = [2,1,3]; %[target,ps_piezo_fast,ps_piezo_slow]
addParameter(p,'ChannelAssignment',defaultChannelAssignment,@isvector);
defaultPeriod = 2; % number of periods expected in one piezo segment. Important for phase computation. 
addParameter(p,'Period',defaultPeriod,@isnumeric);
defaultRecomputeTheta = true;
addParameter(p,'RecomputeTheta',defaultRecomputeTheta,@islogical);
% defaultRecomputeOrth = false;
% addParameter(p,'RecomputeOrth',defaultRecomputeOrth,@islogical);
defaultSavePostselection = true;
addParameter(p,'SavePostselection',defaultSavePostselection,@islogical);
defaultSaveTheta = true;
addParameter(p,'SaveTheta',defaultSaveTheta,@islogical);
defaultSaveOrth = true;
addParameter(p,'SaveOrth',defaultSaveOrth,@islogical);
defaultSelParams = struct('Type','phaseAndAmplitude','Position',[0,0.5,10,1]); % use Type 'phase'?
addParameter(p,'SelectionParameters',defaultSelParams,@isstruct);
parse(p,varargin{:});
c = struct2cell(p.Results);
[chAssign,folder,periodsPerSeg,recomputeTheta,saveps,savetheta,saveOrth,selParams] = c{:};
%[chAssign,folder,periodsPerSeg,recomputeOrth,recomputeTheta,saveOrth,saveps,savetheta,selParams] = c{:};

%% Create folder 'post-data'
if ~exist([pwd 'post-data'],'dir')
    mkdir('post-data')
end

%% Iterate through data files
quantities = struct; % Structure that will contain quantities of interest
tempsaveps = false; % Can change postselection saving behavior per case
selStr = selParamsToStr(selParams);
[filenames,~,positions]= getParametersFromFilenames('Folder',folder,'Parameter','position');
dispstat('','init','timestamp','keepthis',0);

for i = 1:length(filenames)  %1
%     load([folder '\' cell2mat(filenames(fileI))],'X2');
%     [~,n,~] = nPhotons(X2,X2,X2);
    quantities.position(i) = positions(i);
    filename = cell2mat(filenames(i));
    clear X1 X2 X3 theta piezoSign thetaMira
    clear O1 O2 O3 oTheta oThetaMira selX selTheta thetaMiraSel;
    postFilename =  ['post-data/',filename,'-',selStr,'.mat'];
    
    dispstat(['loading quadratures of ' filename],'timestamp','keepthis',0);
    
    if ~saveps && ~recomputeTheta
        try
            load(postFilename);
        catch
            dispstat(['Could not find ',postFilename, ...
                ' loading raw quadratures ...'],'timestamp','keepthis',0);
            load([folder '\' filename]);
            tempsaveps = true;
        end
    else
        load([folder '\' filename]); 
    end
    
    if exist('X1','var')
        % make sure all have the same number of pulses. 
         if ~isequal(size(X1,1),size(X2,1),size(X3,1))
             X1 = X1(1:min([size(X1,1),size(X2,1),size(X3,1)]),:,:);
             X2 = X2(1:min([size(X1,1),size(X2,1),size(X3,1)]),:,:);
             X3 = X3(1:min([size(X1,1),size(X2,1),size(X3,1)]),:,:);
         end
 
        % set which channel ist the target channel etc
        quadratures = zeros([size(X1) 3]);
        quadratures(:,:,:,1) = X1;
        quadratures(:,:,:,2) = X2;
        quadratures(:,:,:,3) = X3;
        Xtg = quadratures(:,:,:,chAssign(1));  
        XpsFast = quadratures(:,:,:,chAssign(2));
        XpsSlow = quadratures(:,:,:,chAssign(3));
        clear('quadratures'); 
    end
    

    %%
    %% Compute Phase and Postselected Variables
    
    [thetaCorr,thetaMiraCorr,XtgCorr,XpsSlowCorr,XpsFastCorr]=corrBoundariesOfPiezosegments(Xtg,XpsSlow,XpsFast,piezoSign);
    %[XtgCorr,XpsSlowCorr,XpsFastCorr]=corrBoundariesOfPiezosegments(Xtg,XpsSlow,XpsFast,piezoSign);
    
    [O1,O2,O3,oTheta,iOrth] = selectOrthogonal(XpsFastCorr,XpsSlowCorr,XtgCorr,thetaCorr,piezoSign);
            oThetaMira = thetaMiraCorr(iOrth);
            

               
        % Compute photon numbers for each channel
        [nTg,nPsFast,nPsSlow] = nPhotons(XtgCorr,XpsFastCorr,XpsSlowCorr); 
        
        %% Postselection
        % We are not sure on which phase we have to postselect. 
        % Here are two possibilites, comment / uncomment one of these possibilites.
        % One possibility is to only postselect on phases in the Husimi function 
        %(which is the phase between diode and LO in the postselection channels):  
                   
%         [selX,selTheta,iSelect] = selectRegion(O1,O2,O3,oTheta,selParams);
%         thetaMiraSel = oThetaMira(iSelect);
        
        %Another possibility is to postselect on a total phase computed
        %by adding the phase of the Husimi function to the phase
        %theta (which is the phase between target channel and postselection
        %channel). Maybe the phase thetaMira, which is the phase in the target channel between LO
        %and Mira from doubleslit should also be added to compute the total
        %phase between diode and Mira from doubleslit in the target
        %channel?
        %Here, we are not sure, whether to add the phases or subtract them.
        % You can try different ways to add phases in
        % selectRegionOfTotalPhase.m

        [selX,selTheta,iSelect] = selectRegionAroundZero(O1,O2,O3,oTheta,oThetaMira,selParams);
        thetaMiraSel = oThetaMira(iSelect);
              

              
        fractionSel = length(selX(:))/length(O1(:));
        quantities.fracSel(i) = fractionSel;
        quantities.lengthSelX(i) = length(selX(:));
        quantities.lengthO1(i) = length(O1(:));
        close all;
        
        
  %  end %if ~exist selX
    
    %% Photon number computation
     % compute photon number of postselected quadratures in the doubleslit.
     % We have to make sure that all phases between LO and Mira signal from
     % doubleslit are included to get a meaningful photon number.
     % Therefore, use uniform sampling of the Mira phases thetaMiraSel. 
    
    [nValues] = deal(zeros(10,1));
    for iN=1:10
        try
            uniformX = seriesUniformSampling(selX,thetaMiraSel,'NBins',100);
        catch
            warning(['Problem with uniformSampling.', ...
                'Use X without uniform sampling.']);
            uniformX = selX;
        end  
        [nValues(iN),~,~] = nPhotons(uniformX,uniformX,uniformX);

    end 
    
    nDs = mean(nValues); % postselected photon numbers in the doubleslit
    quantities.nDs(i) = nDs;
    quantities.nDsStd(i) = std(nValues);

    if (exist('nTg','var'))
        quantities.nTg(i) = nTg;
        quantities.nPsFast(i) = nPsFast;
        quantities.nPsSlow(i) = nPsSlow;
    end

  %% Save workspace variables (because recomputing them takes time)
    % Save Theta
    if savetheta 
            save([folder '\' filename],'thetaCorr','thetaMiraCorr','-append');
    end
 
    if saveOrth
        save([folder '\' filename],'O1','O2','O3','oTheta','oThetaMira','iOrth','-append');
    end
    
    % Save postselected variables
    if saveps || tempsaveps
        save(postFilename, ...
            'selX','selTheta','thetaMiraSel','selParams','nTg','nPsFast','nPsSlow','nDs','fractionSel');      
        tempsaveps = false;
    end
    
    
end

%% Create and write table
save([datestr(date,'yyyy-mm-dd-'),'series3Ch-',selStr,'.mat'],'quantities');
% Load most recent table file 'yyyy-MM-dd-series3Ch.txt'
T = seriesRead3ChTable();
if isempty(T)
    T = cell2table(filenames','VariableNames',{'Filename'});
end

% Update table with new values by looping over 'quantities' variable
fields = fieldnames(quantities);
for iField = 1:numel(fields)
    T.(fields{iField}) = makecol(quantities.(fields{iField}));
end

% write results to a new table file
writetable(T,[datestr(date,'yyyy-mm-dd-'),'series3Ch-',selStr,'.txt']);


%% Plot stuff
plot(quantities.position,quantities.nDs,'o-');
xlabel('Position (mm)');
ylabel('Postselected photon number');
graphicsSettings();
ylim([0.5 2.5]);
savefig([selStr,'-nDsVsPosition.fig']);
print([selStr,'-nDsVsPosition.png'],'-dpng','-r300');
clf();

plot(quantities.position,quantities.nDs./quantities.nTg,'o-');
xlabel('Position (mm)');
ylabel('n_{postselected}/n_{Tg,not postselected} ');
graphicsSettings();
ylim([0 2]);
savefig([selStr,'-nDsNormalizedVsPosition.fig']);
print([selStr,'-nDsNormalizedVsPosition.png'],'-dpng','-r300');
clf();




end 