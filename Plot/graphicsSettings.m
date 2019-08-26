function [] = graphicsSettings
    fontName = 'Times New Roman';
    set(gca,'LineWidth',3,'XColor',[0 0 0], 'YColor', [0 0 0],'Box','on',...
        'FontSize',25,'FontName',fontName,...
        'TickDir','Out');
    set(gca,'DefaultTextInterpreter','latex');
end