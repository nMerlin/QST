function [ dataStruct, Is, nAvs, g2Avs , g2Stds] = g2SeriesFromQuadratures( nResolution, varargin)
%DLSERIES Batch processing of series of g2 measurements, with already
%computed quadratures that are in the folder 'Quadratures'
%Computes g2 either time resolved (set 'G2method', 'time') or photon number
%resolved (set 'G2method', 'bins').
%   'Weight': If set 'yes', g2 from the time resolved method is weighted 
%   with the respective photon numbers. 

% Optional input arguments
%% Validate and parse input arguments
p = inputParser;
defaultG2method = 'time'; % time or bins
addParameter(p,'G2method',defaultG2method);
defaultBins = 100; % to be used with the method 'bins'
addParameter(p,'Bins',defaultBins,@isnumeric);
defaultWeight = 'no'; % to be used with the method 'time'
addParameter(p,'Weight',defaultWeight);
parse(p,varargin{:});
c = struct2cell(p.Results);
[bins,g2method,weight] = c{:};

%% Variables
dataStruct = struct('filename',{},'I',{},'nAv',{},'g2Av',{}, 'g2Std',{});

%% Create data overview
dispstat('','init');
dispstat('Checking filenames ...','timestamp','keepthis');
Contents = dir('Quadratures');
name = {Contents.name};

for iStruct =  1:length(Contents) 
    %get filename
    filename = cell2mat(name(iStruct));
    if strcmp(filename,'.') || strcmp(filename,'..')
        continue
    end
    dataStruct(iStruct).filename = filename;
    % get current or power
    %currentToken = regexpi(filename,'-([0123456789.]*)mA','tokens');
%     currentToken = regexpi(filename,'mW-([0123456789.]*)mW','tokens');
     currentToken = regexpi(filename,'([0123456789,]*)mW-5mW','tokens');
     %currentToken = regexpi(filename,'([0123456789,]*)mW','tokens'); 
     currentToken{1}=strrep(currentToken{1},',','.');
    dataStruct(iStruct).I = str2double(cell2mat(currentToken{1}));
    
    %%% Load data
    cd('Quadratures');
    load(filename);
    
    % Compute g2 values according to g2method
    
    if strcmp(g2method,'time'); 
        [g2vec, ada, times] = g2(X, nResolution);
        nAv = mean(ada);  % Mean photon number  
        if strcmp(weight, 'yes')
            g2Av = g2vec'*ada'/sum(ada);
        else
            g2Av = mean(g2vec);    
        end
        g2Std = sqrt(var(g2vec)); 
        %plot g2
        cd('..');
        mkdir('Plots');
        cd('Plots');
        plotG2Stuff(times, g2vec, ada, [strrep(filename,'.mat','') '-Res-' num2str(nResolution) ]);
        cd('..');
        cd('g2Data');
        save([strrep(filename,'.mat','') '-G2-nResolution-' num2str(nResolution) '-' g2method '.mat'],...
            'times','g2vec','ada');
        cd('..');
    else
        cd('..');
        mkdir('Plots');
        cd('Plots');
        [g2vec, ada, meang2] = g2Bins(X, nResolution, bins, filename); 
        g2Av = meang2;
        g2Std = sqrt(var(g2vec,'omitnan'));  %?? use this or only variance in the middle range?     
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

save(['AverageNandG2-' g2method '.mat'],'Is','nAvs','g2Avs','g2Stds');
xlswrite(['Averages-' g2method '.xls'],[Is' nAvs' g2Avs' g2Stds' ]);
plotNandG2Av(g2method);
end % function

