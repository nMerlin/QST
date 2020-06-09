function [ dataStruct] = plotCrossCorrelationDelaySeries(varargin)
%DLSERIES Pots the crosscorrelations of a 3 Channel series, with already
%computed quadratures that are in the folder 'mat-data'

% Optional input arguments
%% Validate and parse input arguments
p = inputParser;
defaultNSegments = 2;
addParameter(p,'NSegments',defaultNSegments,@isnumeric);
parse(p,varargin{:});
c = struct2cell(p.Results);
[nsegments] = c{:};

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
    % get delay; from the file name format xx-yymm, where xx is the
    % filenumber and yy is the delay which can also be negative and start
    % with a minus sign
    delayToken = regexpi(filename,'([-0123456789,-]*)mm','tokens');
    delay = cell2mat(delayToken{1});
    numberToken = regexpi(delay,'^([0123456789,]*)-','tokens'); 
    number = cell2mat(numberToken{1});
    delay = strrep(delay,[number '-'],'');
    delay = strrep(delay,',','.');
    dataStruct(iStruct).delay = str2double(delay);
    
    %%% Load data
    %cd('Quadratures');
    cd('mat-data');
    load(filename);
    cd('..');
    
    % Plot crosscorrelation
    dispstat('Plot cross corr ...','timestamp','keepthis');
    cd('crosscorrelations');
    [A12,A13,A23] = plotCrossCorrelation(X1, X2, X3, filename, 'Nsegments',nsegments);
    cd('..');
    
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

