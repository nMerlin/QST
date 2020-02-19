function [] = plotNandG2Av(filename,parameter)
    %load('AverageNandG2.mat');
    load(['AverageNandG2-' filename '.mat']);
    plot(min(Is)-2:1:max(Is)+2, ones(length(min(Is)-2:1:max(Is)+2)),'-','linewidth',2,...
        'Color',[1 0.6 0]); %coherent
    hold on;
    plot(min(Is)-2:1:max(Is)+2, 2*ones(length(min(Is)-2:1:max(Is)+2)),'-','linewidth',2,...
        'Color',[132/255 184/255 24/255]); %thermal
    %plot(70*ones(length(0:0.2:3)), 0:0.2:3,'-','linewidth',2,'Color',[0.5 0.5 0.5]); %threshold
    errorbar(Is, g2Avs,g2Stds,'ro','markerSize',7,'markerFaceColor','r','linewidth',1.2); 

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
   
    graphicsSettings;
    grid on; 
    ylim([0.5 2.5]);

    %print('g2overCurrent','-dpng');
    print(['g2overPower-' filename],'-dpng');
    savefig(['g2overPower-' filename '.fig']);

    clf();
    % 
    loglog(Is, nAvs,'o','LineWidth',2);
    ylabel('{<n>}');
    switch parameter
        case 'current'
            xlabel('$I$ (mA)','FontSize',fontsize,'Interpreter','latex');
        case 'power'
            xlabel('Excitation Power (mW)','FontSize',fontsize,'Interpreter','latex');
        case 'delay'
            xlabel('delay','FontSize',fontsize,'Interpreter','latex');
    end
    graphicsSettings;
    grid on; 
    print(['NoverPower-' filename],'-dpng');
    savefig(['NoverPower-' filename '.fig']);
end