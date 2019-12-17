function [] = graphicsSettings
    fontName = 'Times New Roman';
    fontSize = 22;
    set(gca,'LineWidth',3,'XColor',[0 0 0], 'YColor', [0 0 0],'Box','on',...
        'FontSize',fontSize,'FontName',fontName,...
        'TickDir','Out');
    set(gca,'DefaultTextInterpreter','latex');
    grid;
end