function [] = photonNumbersVsPositionWithPostselection(varargin)

%% Validate and parse input arguments
p = inputParser;
defaultFolder = 'mat-data'; % which folder contains the data
addParameter(p,'Folder',defaultFolder);
defaultChannelAssignment = [2,1,3]; %[target,ps_piezo_fast,ps_piezo_slow]
addParameter(p,'ChannelAssignment',defaultChannelAssignment,@isvector);
defaultPeriod = 2; % number of periods expected in one piezo segment. Important for phase computation. 
addParameter(p,'Period',defaultPeriod,@isnumeric);
defaultRecomputeTheta = false;
addParameter(p,'RecomputeTheta',defaultRecomputeTheta,@islogical);
defaultRecomputeOrth = false;
addParameter(p,'RecomputeOrth',defaultRecomputeOrth,@islogical);
defaultSavePostselection = false;
addParameter(p,'SavePostselection',defaultSavePostselection,@islogical);
defaultSaveTheta = false;
addParameter(p,'SaveTheta',defaultSaveTheta,@islogical);
defaultSaveOrth = false;
addParameter(p,'SaveOrth',defaultSaveOrth,@islogical);
defaultSelParams = struct('Type','fullcircle','Position',[2.5,0.5]); % use Type 'phase'?
addParameter(p,'SelectionParameters',defaultSelParams,@isstruct);
parse(p,varargin{:});
c = struct2cell(p.Results);
[chAssign,folder,periodsPerSeg,recomputeOrth,recomputeTheta,saveOrth,saveps,savetheta,selParams] = c{:};

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

for i = 1:length(filenames)
%     load([folder '\' cell2mat(filenames(fileI))],'X2');
%     [~,n,~] = nPhotons(X2,X2,X2);
    quantities.position(i) = positions(i);
    filename = cell2mat(filenames(i));
    clear X1 X2 X3 theta piezoSign
    clear O1 O2 O3 oTheta selX selTheta;
    postFilename =  ['post-data/',filename,'-',selStr,'.mat'];
    
    if ~saveps && ~recomputeTheta
        try
            load(postFilename);
        catch
            dispstat(['Could not find ',postFilename, ...
                ' loading raw quadratures ...'],'timestamp','keepthis',0);
            load([folder '\' filename],'X1','X2','X3','piezoSign');
            tempsaveps = true;
        end
    else
        load([folder '\' filename],'X1','X2','X3','piezoSign');
    end
    
    if exist('X1','var')
        quadratures = zeros([size(X1) 3]);
        quadratures(:,:,:,1) = X1;
        quadratures(:,:,:,2) = X2;
        quadratures(:,:,:,3) = X3;
        Xtg = quadratures(:,:,:,chAssign(1));  % this sets the target channels etc
        XpsFast = quadratures(:,:,:,chAssign(2));
        XpsSlow = quadratures(:,:,:,chAssign(3));
        clear('quadratures'); 
    end
    
    %% Compute Phase and Postselected Variables
    if ~exist('selX','var') % run only if postselection file was not loaded
        if (~exist('theta','var')) || recomputeTheta % run only if theta is not available or should be rewritten
            try
                [theta,~] = computePhase(Xtg,XpsFast,piezoSign,'Period',periodsPerSeg);%brauchen wir die?
            catch
                warning(['Problem using computePhase.', ...
                    ' Assigning random phase.']);
                theta = rand(size(Xtg,1)*size(Xtg,2),size(Xtg,3))*2*pi;
            end
            
        end
        
        %select orthogonal quadratures from the postselection channels
        if (~exist('O1','var')) || recomputeOrth 
            [O1,O2,O3,oTheta,iOrth] = selectOrthogonal(XpsFast,XpsSlow,Xtg,theta,piezoSign);
        end
        
        % Compute photon numbers for each channel
        [nTg,nPsFast,nPsSlow] = nPhotons(Xtg,XpsFast,XpsSlow); 
        
        % postselection from a certain region of the Husimi function                
        [selX,selTheta] = selectRegion(O1,O2,O3,oTheta,selParams);%,'Plot','show','Filename',[filename '-assessTheta']
        fractionSel = length(selX(:))/length(O1(:));
        quantities.fracSel(i) = fractionSel;
        quantities.lengthSelX(i) = length(selX(:));
        quantities.lengthO1(i) = length(O1(:));
        close all;
        
    end %if ~exist selX
    
     % compute photon number of postselected quadratures in the doubleslit
    [nDs,~,~] = nPhotons(selX,selX,selX);
    quantities.nDs(i) = nDs;

    if (exist('nTg','var'))
        quantities.nTg(i) = nTg;
        quantities.nPsFast(i) = nPsFast;
        quantities.nPsSlow(i) = nPsSlow;
    end
    
     %% Save workspace variables (because recomputing them takes time)
    % Save Theta
    if savetheta 
            save([folder '\' filename],'theta','-append');
    end
    
    if saveOrth
        save([folder '\' filename],'O1','O2','O3','oTheta','iOrth','-append');
    end
    
    % Save postselected variables
    if saveps || tempsaveps
        save(postFilename, ...
            'selX','selTheta','selParams','nTg','nPsFast','nPsSlow','nDs','fractionSel');      
        tempsaveps = false;
    end
    
end %file iteration

%% Create and write table
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
save([datestr(date,'yyyy-mm-dd-'),'series3Ch-',selStr,'.mat'],'quantities');

%% Plot stuff
plot(positions,quantities.nDs,'o-');
xlabel('Position (mm)');
ylabel('Postselected photon number');
graphicsSettings();
savefig([selStr,'-nVsPosition.fig']);
print([selStr,'nVsPosition.png'],'-dpng','-r300');
clf();

 %% normalize the interference pattern to absolute maximum and minimum 
% numberOfPeriods = 9;
% y = quantities.nDs;
% [ynew,locs,maxlocs,minlocs] = normalizeInterferencePattern(y,numberOfPeriods);
% plot(positions,ynew,'o');
% graphicsSettings();
% hold on;
% plot(positions,ynew,'b-');
% xlabel('Position (mm)');
% ylabel('Normalized postselected photon number');
% %plot theory curve 
% %% Estimate fit paramters
% y = ynew;
% yMax = mean(y(maxlocs));
% yMin = mean(y(minlocs));
% yRange = (yMax-yMin);
% xphase = -positions(minlocs(1));
% period  = mean(diff(positions(locs)))*2;
% fitFunction = 'a.*(sin(2*pi*x./b + c)).^2';
% [f,~,~] = fit(positions',y,fitFunction, 'StartPoint', [yRange, period*2, xphase] );
% hold on;
% xNew = 0:0.01:positions(end);
% plot(xNew,yRange.*(sin(2*pi*xNew./(f.b) + f.c)).^2 + yMin,'r-','LineWidth',2);
% fittedPeriodInMM = f.b/2;
% savefig([selStr, '-nVsPosition-normalized.fig']);
% print([selStr,'-nVsPosition-normalized.png'],'-dpng','-r300');
% clf();
% 
% %Interferenzwegl?nge in m
% D = 0.6;
% %Wellenl?nge in m
% lambda = 834e-9;
% %Abstand Maxima in m 
% a = fittedPeriodInMM*1e-3;
% %Abstand der Spalte
% dInM = D*lambda/a;


end 