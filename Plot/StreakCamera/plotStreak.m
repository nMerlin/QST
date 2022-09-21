function [decaytime, decaytimeError, Max,Sum] = plotStreak(filenameSIG, filenameBG,varargin)
%PLOTSTREAK plots the picture of a streak camera. It also plots a sum over
%the space over time and makes an exponential decay fit. 
%
%   Input Arguments:
%       filenameSIG: file with the signal data, of '.img' format.
%       It should be located in a folder 'raw-data'.
%       filenameBG: file with background data.
%       'Subtract','yes': subtract background data
%       'Plottype','lin': linear time- and intensity scale for integrated plot
%       'Plottype','log': logarithmic intensity scale for integrated plot

%% Validate and parse input arguments
parser = inputParser;
defaultSubtract = 'yes'; 
addParameter(parser,'Subtract',defaultSubtract);
defaultPlottype = 'lin'; 
addParameter(parser,'Plottype',defaultPlottype);
defaultIntegrationArea = 'full'; 
addParameter(parser,'IntegrationArea',defaultIntegrationArea);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[integrationArea,plottype,subtract] = c{:};

%% load data
    cd('raw-data');
    [M, time] = loadImg(filenameSIG);
    [MBG, ~] = loadImg(filenameBG);
    
    if strcmp(subtract,'yes') % optional: subtract background
        M = M - MBG;
    end
    
%% make surface plot
    x = 0:size(M,2)-1;
    surf(x, time, M);
    colorbar;
    view(0,-90);
    shading flat;
    axis tight;

    fontsize = 24;
    fontName = 'Times New Roman';
    graphicsSettings;
    set(gca,'XGrid','on','YGrid','on');
    xlabel('x (pixel)','FontSize',fontsize,'FontName',fontName);
    ylabel('time (ps)','FontSize',fontsize,'FontName',fontName);
    
    %write overall maximum and sum
    Max = max(max(M));
    Sum = sum(sum(M));
    text(0.1, 0.1, ['max Int ' num2str(Max,'%.0f') ' counts' char(10) ...
        'integr. Int ' num2str(Sum,'%.0f') ' counts'],'FontSize',fontsize-4,...
        'Units','normalized','Color','w');    
    cd('..');
    print([filenameSIG '-subtract-' subtract '-2Dsurfplot.png'],'-dpng','-r300');
    savefig([filenameSIG '-subtract-' subtract '-2Dsurfplot.fig']);
    clf();
    
%% make integrated plot 

     [~,xmaxpixel] = max(max(M));
     line = M(:,xmaxpixel); % pulseshape at the x position with maximum intensity
     [~,timemaxpixel] = max(line);
     cut = M(timemaxpixel,:); % curve shape at the time with maximum intensity
     FWHM = round(fwhm(x,cut)); % spatial width of the curve
     
     switch integrationArea
         case 'full'
            Int = sum(M,2); % summed intensity over complete image
         case 'box'
            Int = sum(M(:,xmaxpixel-FWHM:xmaxpixel+FWHM),2); %summed intensity within 2 times the FWHM   
         case 'line'
            Int = line; % intensity only along a line at the x position with maximum intensity 
     end   
%     if strcmp(subtract,'no')
%         Int = Int - Int(1); %subtract underground; this should be done with a underground measurement
%     end
    
    if strcmp(plottype,'lin')
        plot(time, Int, 'linewidth',2);
        set(gca, 'Ylim',[0 max(Int)+10000]);
        hold on;
        %% fit the falling slope
        [Max, maxIndex] = max(Int);
        timeFit = time(maxIndex:end)- time(maxIndex);
        IntFit = Int(maxIndex:end);
        [f,gof,~] = fit(timeFit,IntFit,'exp1', 'StartPoint',[Max, -0.005]); 
        %f(x) =  a*exp(b*x)
        h = plot(timeFit+time(maxIndex), f.a*exp(f.b*timeFit));
        h.LineWidth = 1.5;
        h.Color = 'r' ; %[0 1 1]
        l = legend('Data','Fit');
        l.FontSize = 22;
        decaytime = -1/f.b;
        level = 2*tcdf(-1,gof.dfe);
        m = confint(f,level); 
        std = m(end,end) - f.b;
        decaytimeError = 1/f.b^2 * std;
        hold off;
        set(gca,'DefaultTextInterpreter','latex');
        text(0.4, 0.5, ['$\tau = $ ' num2str(decaytime,'%.1f') ' $\pm$ ' ...
            num2str(decaytimeError,'%.1f') ' ps'],...
            'FontSize',fontsize-4,...
            'Units','normalized','Interpreter','latex');

    elseif strcmp(plottype,'log')
        semilogy(time(Int>0), Int(Int>0), 'linewidth',2);
        hold on;
        %% fit the falling slope
        [Max, maxIndex] = max(Int);
        IntFit0 = Int(maxIndex:end);
        IntFit = log(IntFit0(IntFit0>0));
        range = max(IntFit)-min(IntFit);
        level = range - 1;
        levelIndex = find(max(IntFit)-IntFit>level,1,'first');
        %IntFit = IntFit(1:levelIndex); %use only values before too much noise
        timeFit = time(maxIndex:end)- time(maxIndex);
        timeFit = timeFit(IntFit0>0);
        %timeFit = timeFit(1:levelIndex);
        [f,gof,~] = fit(timeFit,IntFit,'poly1'); %fity = p1*fitx + p2
        b = f.p1; decaytime = -1/b;
        level = 2*tcdf(-1,gof.dfe);
        m = confint(f,level); 
        std = m(end,1) - b;
        decaytimeError = 1/b^2 * std;
        h = plot(timeFit+time(maxIndex), exp(f.p2)*exp(b*timeFit));
        h.LineWidth = 1.5;
        h.Color = 'r' ; %[0 1 1]
        l = legend('Data','Fit');
        l.FontSize = 22;
        hold off;
        set(gca,'DefaultTextInterpreter','latex');
        text(0.4,0.5, ['$\tau = $ ' num2str(decaytime,'%.1f') ' $\pm$ ' ...
            num2str(decaytimeError,'%.1f') ' ps'],...
            'FontSize',fontsize-4,...
            'Units','normalized','Interpreter','latex');         
    end
    ylabel('Integrated Intensity (a.u.)','FontSize',fontsize,'FontName',fontName);
    xlabel('time (ps)','FontSize',fontsize,'FontName',fontName);
    set(gca,'LineWidth',3,'XColor',[0 0 0], 'YColor', [0 0 0],'Box','on',...
    'FontSize',22,'FontName',fontName,...
    'TickDir','In');
    print([filenameSIG '-subtract-' subtract '-IntegratedPlot-IntArea-' integrationArea '-' plottype '.png'],'-dpng','-r300');
    savefig([filenameSIG '-subtract-' subtract '-IntegratedPlot-IntArea-' integrationArea '-' plottype '.fig']);
    
    %% save the integrated data
    mkdir ('integrated-data');
    cd('integrated-data');
    save([filenameSIG '-subtract-' subtract '-IntArea-' integrationArea '-IntegratedData.mat'],'time','Int');
    cd('..');
    
end