function [ dataStruct, Is, nAvs, g2Avs , g2Stds] = g2SeriesFromQuadratures( nResolution, varargin)
%DLSERIES Batch processing of series of g2 measurements, with already
%computed quadratures that are in the folder 'Quadratures'
%Computes g2 either time resolved (set 'G2method', 'time') or photon number
%resolved (set 'G2method', 'bins').
%   'Weight': If set 'yes', g2 from the time resolved method is weighted 
%   with the respective squared photon numbers. 

% Optional input arguments
%% Validate and parse input arguments
p = inputParser;
defaultG2method = 'time'; % time or bins
addParameter(p,'G2method',defaultG2method);
defaultBins = 100; % to be used with the method 'bins'
addParameter(p,'Bins',defaultBins,@isnumeric);
defaultWeight = 'no'; % to be used with the method 'time'
addParameter(p,'Weight',defaultWeight);
defaultX = 'X'; % which quadrature array from data should be used
addParameter(p,'UseX',defaultX);
defaultParameter = 'power'; % which parameter was changed during the series
addParameter(p,'Parameter',defaultParameter);
parse(p,varargin{:});
c = struct2cell(p.Results);
[bins,g2method,parameter,useX,weight] = c{:};

%% Variables
dataStruct = struct('filename',{},'I',{},'nAv',{},'g2Av',{}, 'g2Std',{});

%% Create data overview
dispstat('','init');
dispstat('Checking filenames ...','timestamp','keepthis');
%Contents = dir('Quadratures');
Contents = dir('mat-data');
name = {Contents.name};

for iStruct =  1:length(Contents) 
    %get filename
    filename = cell2mat(name(iStruct));
    if strcmp(filename,'.') || strcmp(filename,'..') || strcmp(filename,'.txt')
        continue
    end
    dataStruct(iStruct).filename = filename;
    % get current or power
    %currentToken = regexpi(filename,'-([0123456789.]*)mA','tokens');
%     currentToken = regexpi(filename,'mW-([0123456789.]*)mW','tokens');
%      currentToken = regexpi(filename,'([0123456789,]*)mW-5mW','tokens');
    currentToken = regexpi(filename,'([0123456789,]*)mW-4mW','tokens');
     %currentToken = regexpi(filename,'([0123456789,]*)mW','tokens'); 
    % currentToken = regexpi(filename,'mat([0123456789,]*)','tokens');
     currentToken{1}=strrep(currentToken{1},',','.');
    dataStruct(iStruct).I = str2double(cell2mat(currentToken{1}));
    
    if ~exist(['Plots-' useX],'dir')
    mkdir(['Plots-' useX])
    end
    if ~exist(['g2data-' useX],'dir')
    mkdir(['g2data-' useX])
    end
    
    %%% Load data
    %cd('Quadratures');
    cd('mat-data');
    load(filename);
    
    % Compute g2 values according to g2method
    
    switch useX
        case 'X'
        case 'X1'
            X = X1;
        case 'X2'
            X = X2;
        case 'X3'
            X = X3;
    end
    
    if strcmp(g2method,'time'); 
        [g2vec, ada, times] = g2(X, nResolution);
        nAv = mean(ada);  % Mean photon number  
        if strcmp(weight, 'yes')
            g2Av = g2vec'*ada.^2'/sum(ada.^2);
        elseif strcmp(weight, 'slow') 
            g2Av = mean(ada.^2)/mean(ada).^2;
        else
            g2Av = mean(g2vec);    
        end      
        g2Std = sqrt(var(g2vec)); 
        %plot g2
        cd('..');
        cd(['Plots-' useX]);
        plotG2vsTime(times, g2vec, ada, [strrep(filename,'.mat','') '-Res-' num2str(nResolution) ]);
        cd('..');
        cd(['g2Data-' useX]);
        save([strrep(filename,'.mat','') '-G2-nResolution-' num2str(nResolution) '-' g2method '.mat'],...
            'times','g2vec','ada');
        cd('..');
    else
        cd('..');
        cd(['Plots-' useX]);
        [g2vec, ada, meang2, g2Std] = g2Bins(X, nResolution, bins, filename); 
        g2Av = meang2;
        %g2Std = sqrt(var(g2vec,'omitnan'));  %?? use this or only variance in the middle range?     
        nAv = mean(ada,'omitnan'); % Mean photon number  
        cd('..');
        
    end
    dataStruct(iStruct).g2Av = g2Av;
    dataStruct(iStruct).g2Std = g2Std;
    dataStruct(iStruct).nAv = nAv;
    
end % iStruct

Is = cell2mat({dataStruct.I});
nAvs = cell2mat({dataStruct.nAv});
g2Avs = cell2mat({dataStruct.g2Av});
g2Stds = cell2mat({dataStruct.g2Std});

save(['AverageNandG2-' useX '-' g2method '-weight-' weight '.mat'],'Is','nAvs','g2Avs','g2Stds');
xlswrite(['Averages-' useX '-' g2method '-weight-' weight '.xls'],[Is' nAvs' g2Avs' g2Stds' ]);
plotNandG2Av([useX '-' g2method '-weight-' weight],parameter);
end % function

