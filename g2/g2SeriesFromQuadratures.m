function [ dataStruct, Is, nAvs, g2Avs , g2Stds] = g2SeriesFromQuadratures( nResolution)
%DLSERIES Batch processing of diode laser current series, with already
%computed quadratures that are in the folder 'Quadratures'
%   Detailed explanation goes here

% Optional input arguments


%%% Variables
dataStruct = struct('filename',{},'I',{},'nAv',{},'g2Av',{}, 'g2Std',{});

%%% Create data overview
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
     currentToken{1}=strrep(currentToken{1},',','.');
    dataStruct(iStruct).I = str2double(cell2mat(currentToken{1}));
    
    %%% Load data
    cd('Quadratures');
    load(filename);
    
    % Compute g2 values
    [g2vec, ada, times] = g2(X, nResolution);
    g2Av = mean(g2vec);
    g2Std = sqrt(var(g2vec));
    dataStruct(iStruct).g2Av = g2Av;
    dataStruct(iStruct).g2Std = g2Std;
    % Mean photon number  
    nAv = mean(ada);
    dataStruct(iStruct).nAv = nAv;
    
    %plot g2
    cd('..');
    cd('Plots');
    plotG2Stuff(times, g2vec, ada, [strrep(filename,'.mat','') '-Res-' num2str(nResolution) ]);
    cd('..');
    cd('g2Data');
    save([strrep(filename,'.mat','') '-G2-nResolution-' num2str(nResolution) '.mat'],...
        'times','g2vec','ada');
    cd('..');
    
end % iStruct

Is = cell2mat({dataStruct.I});
nAvs = cell2mat({dataStruct.nAv});
g2Avs = cell2mat({dataStruct.g2Av});
g2Stds = cell2mat({dataStruct.g2Std});

save('AverageNandG2.mat','Is','nAvs','g2Avs','g2Stds');
xlswrite('Averages.xls',[Is' nAvs' g2Avs' g2Stds' ]);
end % function

