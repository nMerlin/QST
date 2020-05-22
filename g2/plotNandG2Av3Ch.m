function [] = plotNandG2Av3Ch(filename,parameter)
    %load('AverageNandG2.mat');
    load(['AverageNandG2-X1-' filename '.mat']);
    Is1 = Is;
    nAvs1 = nAvs;
    g2Avs1 = g2Avs;
    g2Stds1 = g2Stds;
    load(['AverageNandG2-X2-' filename '.mat']);
    Is2 = Is;
    nAvs2 = nAvs;
    g2Avs2 = g2Avs;
    g2Stds2 = g2Stds;
    load(['AverageNandG2-X3-' filename '.mat']);
    Is3 = Is;
    nAvs3 = nAvs;
    g2Avs3 = g2Avs;
    g2Stds3 = g2Stds;
    
    errorbar(Is1, g2Avs1,g2Stds1,'ro','markerSize',7,'markerFaceColor','r','markerEdgeColor','r','linewidth',1.2); 
    hold on;
    errorbar(Is2, g2Avs2,g2Stds2,'ko','markerSize',7,'markerFaceColor','k','markerEdgeColor','k','linewidth',1.2);
    errorbar(Is3, g2Avs3,g2Stds3,'bo','markerSize',7,'markerFaceColor','b','markerEdgeColor','b','linewidth',1.2);
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
    legend('Channel 1','Channel 2','Channel 3','location','northwest');
    graphicsSettings;
    ylim([0.5 2.5]);

    %print('g2overCurrent','-dpng');
    print(['g2overPower-3Ch-' filename],'-dpng');
    savefig(['g2overPower-3Ch-' filename '.fig']);

    clf();
    % 
    loglog(Is1, nAvs1,'o','LineWidth',2);
    hold on;
    loglog(Is2, nAvs2,'o','LineWidth',2);
    loglog(Is3, nAvs3,'o','LineWidth',2);
    ylabel('{<n>}');
    switch parameter
        case 'current'
            xlabel('$I$ (mA)','FontSize',fontsize,'Interpreter','latex');
        case 'power'
            xlabel('Excitation Power (mW)','FontSize',fontsize,'Interpreter','latex');
        case 'delay'
            xlabel('delay','FontSize',fontsize,'Interpreter','latex');
    end
    legend('Channel 1','Channel 2','Channel 3','location','northwest');
    graphicsSettings
    print(['NoverPower-3Ch-' filename],'-dpng');
    savefig(['NoverPower-3Ch-' filename '.fig']);
end