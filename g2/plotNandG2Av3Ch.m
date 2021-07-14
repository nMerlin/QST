function [] = plotNandG2Av3Ch(filename,parameter, varargin)
%% Validate and parse input arguments
p = inputParser;
defaultChNumber = 2; % How many channels should be plotted 
addParameter(p,'ChNumber',defaultChNumber,@isnumeric);
parse(p,varargin{:});
c = struct2cell(p.Results);
[ChNumber] = c{:};
%load('AverageNandG2.mat');
    load(['AverageNandG2-X1-' filename '.mat']);
    Is1 = Is;
    nAvs1 = nAvs;
    g2Avs1 = g2Avs;
    g2Stds1 = g2Stds;
    if ChNumber == 2|| ChNumber ==3
        load(['AverageNandG2-X2-' filename '.mat']);
        Is2 = Is;
        nAvs2 = nAvs;
        g2Avs2 = g2Avs;
        g2Stds2 = g2Stds;
    end
    if ChNumber == 3
        load(['AverageNandG2-X3-' filename '.mat']);
        Is3 = Is;
        nAvs3 = nAvs;
        g2Avs3 = g2Avs;
        g2Stds3 = g2Stds;
    end
    
    l = errorbar(Is1, g2Avs1,g2Stds1,'ko','markerSize',7,'markerFaceColor','k','markerEdgeColor','k','linewidth',1.2);
    f = gca;
    f.XScale = 'log';
    if ChNumber ==2 || ChNumber ==3
        l.DisplayName='Channel 1';
        hold on;
        errorbar(Is2, g2Avs2,g2Stds2,'ro','markerSize',7,'markerFaceColor','r','markerEdgeColor','r','linewidth',1.2,'DisplayName','Channel 2');
    end
    if ChNumber == 3    
        errorbar(Is3, g2Avs3,g2Stds3,'bo','markerSize',7,'markerFaceColor','b','markerEdgeColor','b','linewidth',1.2,'DisplayName','Channel 3');
    end
    fontsize =22;
    ylabel('$ \langle g^{(2)}(0) \rangle $','FontSize',fontsize, 'Interpreter','latex');
    switch parameter
        case 'current'
            xlabel('$I$ (mA)','FontSize',fontsize,'Interpreter','latex');
        case 'power'
            xlabel('Excitation Power (mW)','FontSize',fontsize,'Interpreter','latex');
        case 'delay'
            xlabel('delay','FontSize',fontsize,'Interpreter','latex');
    end
    legend('location','best');
    graphicsSettings;
    ylim([0.5 2.5]);

    %print('g2overCurrent','-dpng');
    print(['g2overPower-allCh-ChNumber-' num2str(ChNumber) filename],'-dpng');
    savefig(['g2overPower-allCh-ChNumber-' num2str(ChNumber) filename '.fig']);

    clf();
    % 
    loglog(Is1, nAvs1,'ko','LineWidth',2,'DisplayName','Channel 1');
    hold on;
    loglog(Is2, nAvs2,'ro','LineWidth',2,'DisplayName','Channel 2');
    if ChNumber == 3 
        loglog(Is3, nAvs3,'bo','LineWidth',2,'DisplayName','Channel 3');
    end
    ylabel('{<n>}');
    switch parameter
        case 'current'
            xlabel('$I$ (mA)','FontSize',fontsize,'Interpreter','latex');
        case 'power'
            xlabel('Excitation Power (mW)','FontSize',fontsize,'Interpreter','latex');
        case 'delay'
            xlabel('delay','FontSize',fontsize,'Interpreter','latex');
    end
    legend('location','northwest');
    graphicsSettings
    print(['NoverPower-allCh-ChNumber-' num2str(ChNumber) filename],'-dpng');
    savefig(['NoverPower-allCh-ChNumber-' num2str(ChNumber) filename '.fig']);
end