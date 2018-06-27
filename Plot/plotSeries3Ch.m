function plotSeries3Ch(T,varargin)
%PLOTSERIES3CH Plot results from table T created by SERIES3CH.

%% Validate and parse input arguments
p = inputParser;
defaultFilename = '';
addParameter(p,'Filename',defaultFilename,@isstr);
defaultType = 'Delay';
addParameter(p,'Type',defaultType,@isstr);
parse(p,varargin{:});
c = struct2cell(p.Results);
[filename,typestr] = c{:};

fig = figure;
hold on;
switch typestr
    case 'DelayMeanVarX'
        %% Create plot delays versus postselected variance
        % Plot data
        xAxis = T.Delay;
        varX = T.discMeanVar;
        plot(xAxis,varX,'o','DisplayName','(\Delta Q)^2');
        ax = get(fig,'CurrentAxes');
        xlabel(ax,'Delay (fs)');
        
        % Fit with Gaussian
        fo = fitoptions('Method','NonlinearLeastSquares', ...
            'StartPoint',[-5 0 600 8.3]);
        ft = fittype('a*exp(-(x-b)^2/(2*s^2))+c','options',fo);
        res = fit(xAxis,varX,ft);
        plot(ax,min(xAxis):1:max(xAxis),res(min(xAxis):1:max(xAxis)),'r', ...
            'DisplayName',['Gaussian with \sigma=', ...
            num2str(round(res.s)),' fs']);
        xlim(ax,[min(xAxis) max(xAxis)]);
        
        % Lower Bound
        %minVar = compute3ChLimit(T.nX2,T.nX3,T.nX1);
        %plot(xAxis,minVar,'x','DisplayName',['Limit ca. ', ...
        %    num2str(min(minVar))]);

        % Upper Bound
        %plot(xAxis,T.nX1+0.5,'*','DisplayName','Upper Bound');

        % Labels
        ylabel(ax,'Postselected Variance');
    case 'DelayDiscAmpl'
        %% Create plot delays versus amplitude (discretization method)
        xAxis = T.Delay;
        yAxis = T.discAmpl;
        plot(xAxis,yAxis,'o','DisplayName','|\alpha|');
        ax = get(fig,'CurrentAxes');
        xlabel(ax,'Delay (fs)');
        ylabel(ax,'Coherent Amplitude');
        xlim(ax,[min(xAxis) max(xAxis)]);
        
        % Fit with Gaussian
        fo = fitoptions('Method','NonlinearLeastSquares', ...
            'StartPoint',[-5 0 600 8.3]);
        ft = fittype('a*exp(-(x-b)^2/(2*s^2))+c','options',fo);
        res = fit(xAxis,yAxis,ft);
        fitx = min(xAxis):1:max(xAxis);
        fity = res(min(xAxis):1:max(xAxis));
        plot(ax,fitx,fity,'r','DisplayName',['Gaussian; \sigma=', ...
            num2str(round(abs(res.s))),' fs']);
    case 'Photons'
        %% Create plot for series with different photon numbers
        nX1 = T.nX1; nX2 = T.nX2; nX3 = T.nX3;
        xAxis = nX1./(nX2+nX3);
        varX = T.discMeanVar;
        plot(xAxis,varX,'o','DisplayName', ...
            ['Minimum: ',num2str(min(varX))]);
        ax = get(fig,'CurrentAxes');
        xlabel(ax,'n_t/n_{ps}');
end

%% Common figure manipulation
set(fig,'Color','w');
legend(ax,'show');
hold off;

%% Write figure to file and close it
if ~isempty(filename)
    savefig(fig,filename);
    close all;
end

end

