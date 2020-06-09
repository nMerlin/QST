function [] = graphicsSettings
    fontName = 'Times New Roman';
    fontSize = 22;
    h = findobj(gca,'Type','line');
    set(h,'linewidth',2);
    set(gca,'linewidth',3,'XColor',[0 0 0], 'YColor', [0 0 0],'Box','on',...
        'FontSize',fontSize,'FontName',fontName,...
        'TickDir','Out');
    set(gca,'DefaultTextInterpreter','latex');
    grid;
    dcmObj = datacursormode;  % Turn on data cursors and return the
                          %   data cursor mode object
    set(dcmObj, 'UpdateFcn', @dataTipUpdateFcn); % increases the displayed digits of data tips 
end