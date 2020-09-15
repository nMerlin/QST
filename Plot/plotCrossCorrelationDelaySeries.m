function [ dataStruct] = plotCrossCorrelationDelaySeries(varargin)
%DLSERIES Pots the crosscorrelations of a 3 Channel series, with already
%computed quadratures that are in the folder 'mat-data'

% Optional input arguments
%% Validate and parse input arguments
p = inputParser;
defaultNSegments = 10;
addParameter(p,'NSegments',defaultNSegments,@isnumeric);
defaultParameter = 'delay'; % which parameter was changed during the series
addParameter(p,'Parameter',defaultParameter);
parse(p,varargin{:});
c = struct2cell(p.Results);
[nsegments,parameter] = c{:};

%% Variables
dataStruct = struct('filename',{},'delay',{},'A12',{},'A13',{}, 'A23',{});

%% Create data overview
dispstat('','init');
dispstat('Checking filenames ...','timestamp','keepthis');
%Contents = dir('Quadratures');
Contents = dir('mat-data');
name = {Contents.name};

if ~exist('crosscorrelations','dir')
    mkdir('crosscorrelations')
end

for iStruct =  1:length(Contents) 
    %get filename
    filename = cell2mat(name(iStruct));
    if strcmp(filename,'.') || strcmp(filename,'..') || strcmp(filename,'.txt')
        continue
    end
    dataStruct(iStruct).filename = filename;
    
    switch parameter
        case 'power'    
            %currentToken = regexpi(filename,'([0123456789,]*)mW','tokens');
            currentToken = regexpi(filename,'([0123456789,]*)mW-4mW','tokens');
             currentToken{1}=strrep(currentToken{1},',','.');
             dataStruct(iStruct).delay = str2double(cell2mat(currentToken{1}));
        case 'delay'
            delayToken = regexpi(filename,'([-0123456789,-]*)mm','tokens');
            delay = cell2mat(delayToken{1});
            numberToken = regexpi(delay,'^([0123456789,]*)-','tokens'); 
            number = cell2mat(numberToken{1});
            delay = strrep(delay,[number '-'],'');
            delay = strrep(delay,',','.');
            delay = str2double(delay);
            c = 299792458; % in m/s
            delay = 2*delay/1000/c*10^12; %delay in ps   
            dataStruct(iStruct).delay = delay;
        case 'no' 
            dataStruct(iStruct).delay = 0;
    end
   
    
    %%% Load data
    %cd('Quadratures');
    cd('mat-data');
    load(filename);
    cd('..');
    
    % Plot crosscorrelation
    dispstat('Plot cross corr ...','timestamp','keepthis');
    [A12,A13,A23] = plotCrossCorrelation(X1, X2, X3, ['crosscorrelations\' filename], ...
        'Nsegments',nsegments);
    
    dataStruct(iStruct).A12 = A12;
    dataStruct(iStruct).A13 = A13;
    dataStruct(iStruct).A23 = A23;

end % iStruct
delay = cell2mat({dataStruct.delay});
A12 = cell2mat({dataStruct.A12});
A13 = cell2mat({dataStruct.A13});
A23 = cell2mat({dataStruct.A23});

c = 299792458; % in m/s
time = 2*delay/1000/c; %delay in mm

cd('crosscorrelations');
save('CrossCorrDelaySeries.mat','delay','time','A12','A13','A23');

plot(delay,A12,'o-','linewidth',3); hold on;
plot(delay,A13,'o-','linewidth',3);
plot(delay,A23,'o-','linewidth',3);
legend('X1*X2','X1*X3','X2*X3');
graphicsSettings;
xlabel('delay (mm)');
ylabel('Amplitude');
savefig('Cross-Correlations-vs-Delay.fig');
print('Cross-Correlations-vs-Delay.png','-dpng','-r300');
clf();

plot(time*10^12,A12,'o-','linewidth',3); hold on;
plot(time*10^12,A13,'o-','linewidth',3);
plot(time*10^12,A23,'o-','linewidth',3);
legend('X1*X2','X1*X3','X2*X3');
graphicsSettings;
xlabel('time delay (ps)');
ylabel('visibility');
savefig('Cross-Correlations-vs-time.fig');
print('Cross-Correlations-vs-time.png','-dpng','-r300');
clf();
cd('..');

end % function

