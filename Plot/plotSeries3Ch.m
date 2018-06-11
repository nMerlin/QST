function plotSeries3Ch(T,varargin)
%PLOTSERIES3CH Plot results from SERIES3CH table T.

%% Validate and parse input arguments
p = inputParser;
defaultSave = false;
addParameter(p,'Save',defaultSave,@islogical);
defaultType = 'Delay';
addParameter(p,'Type',defaultType,@isstr);
parse(p,varargin{:});
c = struct2cell(p.Results);
[saveflag,typestr] = c{:};

switch typestr
    case 'Delay'
        %% Create plot for delay line series
        xAxis = T.Delay;
        varX = T.meanVarX;
        plot(xAxis,varX,'o','DisplayName', ...
            ['Minimum: ',num2str(min(varX))]);
        xlabel('Delay (fs)');
        
        % Fit with Gaussian
        fo = fitoptions('Method','NonlinearLeastSquares', ...
            'StartPoint',[-5 0 600 8.3]);
        ft = fittype('a*exp(-(x-b)^2/t^2)+c','options',fo);
        res = fit(xAxis,varX,ft);
        hold on;
        plot(min(xAxis):1:max(xAxis),res(min(xAxis):1:max(xAxis)),'r', ...
            'DisplayName','Fit');
        xlim([min(xAxis) max(xAxis)]);
        hold off;
    case 'Photons'
        %% Create plot for series with different photon numbers
        nX1 = T.nX1; nX2 = T.nX2; nX3 = T.nX3;
        xAxis = nX1./(nX2+nX3);
        varX = T.meanVarX;
        plot(xAxis,varX,'o','DisplayName', ...
            ['Minimum: ',num2str(min(varX))]);
        xlabel('n_t/n_{ps}');
end

%% Plot limits & labels
hold on;
% Lower Bound
minVar = compute3ChLimit(T.nX2,T.nX3,T.nX1);
plot(xAxis,minVar,'x','DisplayName',['Limit ca. ',num2str(min(minVar))]);

% Upper Bound
plot(xAxis,T.nX1+0.5,'*','DisplayName','Upper Bound');

% Labels
ylabel('Postselected Variance');
legend('show');
title('3-Channel Series');
hold off;

%% Write figure to file
if saveflag
    saveA5Landscape([datestr(date,'yyyy-mm-dd-'),'plotSeries3Ch']);
end

end

