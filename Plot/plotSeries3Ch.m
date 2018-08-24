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
formatFigA5(fig);
legendLocation = 'northeast';
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
        
        % Add fit results
        fitString = evalc('disp(res)');
        text(min(xlim),mean(ylim),fitString);
        
        % Lower Bound
        %minVar = compute3ChLimit(T.nX2,T.nX3,T.nX1);
        %plot(xAxis,minVar,'x','DisplayName',['Limit ca. ', ...
        %    num2str(min(minVar))]);

        % Upper Bound
        %plot(xAxis,T.nX1+0.5,'*','DisplayName','Upper Bound');

        % Labels
        ylabel(ax,'Postselected Variance');
        legendLocation = 'southeast';
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
        
        % Add fit results
        fitString = evalc('disp(res)');
        text(min(xlim),mean(ylim),fitString);
        
    case 'Photons'
        %% Create plot for series with different photon numbers
        nX1 = T.nX1; nX2 = T.nX2; nX3 = T.nX3;
        nSum = nX1 + nX2 + nX3;
        xAxis = nX1./(nX1+nX2+nX3);
        varX = T.discMeanVar;
        xTheo = linspace(0,1,100);
        minVar = (1/mean(nSum)+1+xTheo)./(2*(1/mean(nSum)+1-xTheo));
        maxVar = xTheo*mean(nSum)+0.5;
        hold on;
        plot(xTheo,minVar,'DisplayName','Theory Phase-Sensitive', ...
            'LineWidth',2);
        plot(xTheo,maxVar,'DisplayName','Theory Phase-Randomized', ...
            'LineWidth',2);
        plot(xAxis,varX,'o','DisplayName','Measured','MarkerSize',10, ...
            'MarkerEdgeColor','k','MarkerFaceColor','w','LineWidth',2);
        hold off;
        ax = get(fig,'CurrentAxes');
        xlabel('n_t/(n_t+n_{ps})','FontSize',26);
        ylabel('Var(Q_t)','FontSize',26);
        box on;
        legendLocation = 'NorthWest';
        ax.FontSize = 22;
end

%% Common figure manipulation
set(fig,'Color','w');
legend(ax,'Location',legendLocation);
legend(ax,'show');
if ~isempty(filename)
    selParams = selStrToParams(filename);
    if ~isempty(selParams)
        title(ax,selParamsToStr(selParams));
    end
end
hold off;

%% Write figure to file and close it
if ~isempty(filename)
    savefig(fig,filename);
    close all;
end

end

