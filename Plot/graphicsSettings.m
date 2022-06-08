function [] = graphicsSettings(varargin)
p = inputParser;
defaultFontsize = 19;
addParameter(p,'Fontsize',defaultFontsize,@isnumeric);
parse(p,varargin{:});
c = struct2cell(p.Results);
[fontSize] = c{:};


%     fontName = 'Times New Roman';
    fontName = 'Arial';
    h = findobj(gca,'Type','line');
    ax = gca;
    set(h,'linewidth',2);
    set(ax,'linewidth',2,'XColor',[0 0 0], 'YColor', [0 0 0],'Box','on',...
        'FontSize',fontSize,'FontName',fontName,...
        'TickDir','in');
    %set(ax,'DefaultTextInterpreter','latex');
    grid;
    dcmObj = datacursormode;  % Turn on data cursors and return the
                          %   data cursor mode object
    set(dcmObj, 'UpdateFcn', @dataTipUpdateFcn); % increases the displayed digits of data tips 
  
    %change exponent of axis
%     y = ax.YAxis;
%     set(y,'Exponent',4);
    
    %reduce white margins around figure 
%     outerpos = ax.OuterPosition;
%     ti = ax.TightInset;
%     left = outerpos(1) + ti(1);
%     bottom = outerpos(2) + ti(2);
%     ax_width = outerpos(3)- ti(1) - ti(3);
%     ax_height = outerpos(4) - ti(2) - ti(4);
%     set(gca,'Position',[left bottom ax_width ax_height]);
%     fig = gcf;
%     fig.PaperPositionMode = 'auto';
%     fig_pos = fig.PaperPosition;
%     fig.PaperSize = [fig_pos(3) fig_pos(4)];

end