function [correlations, delayArray] = statisticalIndependence(filename, varargin) 
%This function computes and plots the correlation between data from
%FILENAME and the data shifted about a delay that is varied. 
%PULSES is the maximum number of pulses about which the data should be
%shifted.
%TYPE is the kind of data you want to use. You can either choose the
%raw data with 'rawdata' or the quadratures (integration of raw data) with
% 'quadratures'. When choosing rawdata, the computation may take long.

%Output arguments:
%CORRELATIONS contains the correlation for all the delays and all three
%channels.
%DELAYARRAY contains the delays.

%% Validate and parse input arguments
p = inputParser;
defaultPulses = 20;  %1.2
% Choose 'plot' for a graphical output
defaultPlot = 'plot';
defaultType = 'quadratures';
addParameter(p,'Pulses',defaultPulses,@isnumeric);
addParameter(p,'Plot',defaultPlot);
addParameter(p,'Type',defaultType);
parse(p,varargin{:});
c = struct2cell(p.Results);
[plotArg, pulses,type] = c{:};

%% load data
[data8bit,config,~]=load8BitBinary(filename,'dontsave');

%% get data to be used
 
if strcmp(type,'rawdata')  %use rawdata
    X1 = data8bit(:,:,1);
    X2 = data8bit(:,:,2);
    X3 = data8bit(:,:,3);
    [segmentlength, segments] = size(X1);
    X1 = single(X1(1:round(segmentlength/2),1:round(segments/4)));
    X2 = single(X2(1:round(segmentlength/2),1:round(segments/4)));
    X3 = single(X3(1:round(segmentlength/2),1:round(segments/4)));
    maxDelay = round(segmentlength/1000*pulses);
    delayArray = 1:2:maxDelay;
end

if strcmp(type,'quadratures')  %use quadratures
    X = computeQuadratures(data8bit,config,1);
    X1 = X(:,:,1);
    X2 = X(:,:,2);
    X3 = X(:,:,3);
    maxDelay = pulses;
    delayArray = 1:maxDelay;
end

%% compute correlations
correlations = zeros(3,length(delayArray));
tic;
for i = 1:length(delayArray)
    correlations(1,i) = delayCorr(X1,delayArray(i));
    correlations(2,i) = delayCorr(X2,delayArray(i));
    correlations(3,i) = delayCorr(X3,delayArray(i));
end
toc;

%% plot
if strcmp(plotArg,'plot')
    plot(delayArray, correlations(1,:),'o','MarkerEdgeColor','b','MarkerFaceColor','b',...
        'MarkerSize',13);
    hold on;
    plot(delayArray, correlations(2,:),'d','MarkerEdgeColor','r','MarkerFaceColor','r',...
        'MarkerSize',13);
    plot(delayArray, correlations(3,:),'s','MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerSize',13);
    legend('Channel 1', 'Channel 2','Channel 3','location','best');
    ylabel('correlation coefficient');
    
    if strcmp(type,'quadratures')
        title('Correlation Coefficient of Pulses');
        xlabel('Pulse distance');
    end
    
    if strcmp(type,'rawdata')
        title('Correlation Coefficient of Samples');
        xlabel('Sample distance');
    end
   
end

end


function corr = delayCorr(X, delay)
    A = X(1:end-delay);
    B = X(1+delay:end);
    corr = corrcoef(A,B);
    corr = corr(1,2);
end
