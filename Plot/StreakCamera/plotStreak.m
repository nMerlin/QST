function [decaytime, decaytimeError] = plotStreak(filenameSIG, filenameBG,varargin)
%PLOTSTREAK plots the picture of a streak camera. It also plots a sum over
%the space over time and makes an exponential decay fit. 
%
%   Input Arguments:
%       filenameSIG: file with the signal data, of '.dat' format.
%       It should be located in a folder 'raw-data'.
%       filenameBG: file with background data.

%% Validate and parse input arguments
parser = inputParser;
defaultSubtract = 'yes'; % Rotation angle in degrees
addParameter(parser,'Subtract',defaultSubtract);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[subtract] = c{:};

%% load data
    %datafile = '01-signal-Gain0-Exposure100.dat';
    cd('raw-data');
    M=load(filenameSIG);
    MBG = load(filenameBG);
    
    if strcmp(subtract,'yes') % optional: subtract background
        M = M - MBG;
    end
    
%% make surface plot
    %time = 0:size(M,1)-1;
    prf = dlmread('profil.prf',',',5,0);
    time = prf(:,1);
    x = 0:size(M,2)-1;
    surf(x, time, M);
    colorbar;
    view(0,-90);
    shading flat;
    axis tight;

    fontsize = 24;
    fontName = 'Times New Roman';
    set(gca,'LineWidth',3,'XColor',[0 0 0], 'YColor', [0 0 0],'Box','on',...
        'FontSize',22,'FontName',fontName,...
        'TickDir','Out');
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
    print([filenameSIG '-2Dsurfplot.png'],'-dpng','-r300');
    clf();
    
%% make integrated plot 
    Int = sum(M,2);
    if strcmp(subtract,'no')
        Int = Int - Int(end); %subtract underground; this should be done with a underground measurement
    end
    plot(time, Int, 'linewidth',2);
    ylabel('Integrated Intensity (a.u.)','FontSize',fontsize,'FontName',fontName);
    xlabel('time (ps)','FontSize',fontsize,'FontName',fontName);
    set(gca,'LineWidth',3,'XColor',[0 0 0], 'YColor', [0 0 0],'Box','on',...
        'FontSize',22,'FontName',fontName,...
        'TickDir','In', 'Ylim',[0 max(Int)+10000]);
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
    print([filenameSIG '-IntegratedPlot.png'],'-dpng','-r300');
end