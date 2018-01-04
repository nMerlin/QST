function [correlations, delayArray] = statisticalIndependence(data, varargin) 
% This function computes and plots the autocorrelation of the given data.
%
% Input Arguments:
%   DATA: 1D,2D or 3D matrix; the autocorrelation will be calculated along
%       the first dimension.
%
% Output Arguments:
%   CORRELATIONS: Contains the correlation for all the delays and all three
%       channels.
%   DELAYARRAY: Is a vector with all employed delays.
%
% Optional Input Arguments ('Name',Value):
%   'Pulses',N: The autocorrelation will be computed for the shifts
%       1,2,...,N. Default: N=20.
%   'Plot','Val': A plot of the resulting correlation values can be shown,
%       if 'Val'='show'. Default: 'Val'='show'.
%   'Type','Val': Specifies the type of the input data matrix. Default is
%       'Val'='preformatted', not modifying the DATA matrix. When using a
%       raw data DATA-matrix, employ the option 'rawdata'.
%   'Method','Val': 

%% Validate and parse input arguments
p = inputParser;
defaultPulses = 20;
defaultPlot = 'show';
defaultType = 'preformatted';
defaultMethod = 'matrix';
addParameter(p,'Pulses',defaultPulses,@isnumeric);
addParameter(p,'Plot',defaultPlot);
addParameter(p,'Type',defaultType);
addParameter(p,'Method',defaultMethod);
parse(p,varargin{:});
c = struct2cell(p.Results);
[method, plotArg, pulses,type] = c{:};

%% get data to be used
 
if strcmp(type,'rawdata')  %use rawdata
    X1 = data(:,:,1);
    X2 = data(:,:,2);
    X3 = data(:,:,3);
    [segmentlength, segments] = size(X1);
    X1 = single(X1(1:round(segmentlength/2),1:round(segments/4)));
    X2 = single(X2(1:round(segmentlength/2),1:round(segments/4)));
    X3 = single(X3(1:round(segmentlength/2),1:round(segments/4)));
    maxDelay = round(segmentlength/1000*pulses);
    delayArray = 1:maxDelay;
    [locs,~] = pointwiseVariance(X1);
    start = locs(1);
end

if strcmp(type,'preformatted')  %use quadratures
    % Workaround for single-channel data arrays:
    if size(data,3)==1
        data = repmat(data,1,1,3);
    end
    X1 = data(:,:,1);
    X2 = data(:,:,2);
    X3 = data(:,:,3);
    maxDelay = pulses;
    delayArray = 1:maxDelay;
    start = 1;
end

%% compute correlations
correlations = zeros(3,length(delayArray));
tic;
for i = 1:length(delayArray)
    
    if strcmp(method,'matrix')  %use matrixwise computation of correlation
        correlations(1,i) = delayCorr(X1,delayArray(i));
        correlations(2,i) = delayCorr(X2,delayArray(i));
        correlations(3,i) = delayCorr(X3,delayArray(i));
    end
    
    if strcmp(method,'vector')
        correlations(1,i) = delayCorrVectorwise(X1,delayArray(i),start);
        correlations(2,i) = delayCorrVectorwise(X2,delayArray(i),start);
        correlations(3,i) = delayCorrVectorwise(X3,delayArray(i),start);
    end
    
end
toc;

%% plot
if strcmp(plotArg,'show')
    plot(delayArray, correlations(1,:),'o','MarkerEdgeColor','b','MarkerFaceColor','b',...
        'MarkerSize',13);
    hold on;
    plot(delayArray, correlations(2,:),'d','MarkerEdgeColor','r','MarkerFaceColor','r',...
        'MarkerSize',13);
    plot(delayArray, correlations(3,:),'s','MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerSize',13);
    legend('Channel 1', 'Channel 2','Channel 3','location','best');
    ylabel('correlation coefficient');
    if strcmp(type,'preformatted')
        title('Correlation Coefficient of Pulses');
        xlabel('Pulse distance');
    end
    
    if strcmp(type,'rawdata')
        plot(locs(1:pulses)-start,0*ones(pulses),'v','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',9);
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

function corr = delayCorrVectorwise(X, delay, start)
    A = X(start,:);
    B = X(start+delay,:);
    corr = corrcoef(A,B);
    corr = corr(1,2);

%     Averaging over several vectors.
%     corrVector = zeros(length(start:size(X,1)-delay),1);
%     for s = start:size(X,1)-delay
%         A = X(s,:);
%         B = X(s+delay,:);
%         corrMatrix = corrcoef(A,B);
%         corrVector(s-start+1) = corrMatrix(1,2);
%     end
%     corr = mean(corrVector);
end
