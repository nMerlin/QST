function plotWigner(WF,varargin)
%PLOTWIGNER Plots given Wigner Function
%
% Optional Input Arguments:
%   'Filename': Exports the figure to a pdf with the given filename.
%       Default is '' and saves nothing.
%   'Handle': With the default [] it creates a new plot. Axis handle to
%       plot into. Not implemented for 'Image'.
%   'PQ': Specify p axis. Default is -20:0.125:20.
%   'Style': The 'advanced' style provides a 3D plot with 2D projection
%       beneath and 2D projections on the q and p axes. The '2D' style
%       provides a 2D projection. Default is '3D' with 3D surface plot.
%   'ZLim': Limits of z-axis. Default is [].
%
% Notes:
%   p and q have to be set manually

%% Validate and parse input arguments
p = inputParser;
defaultCBar = 'off';
addParameter(p,'ColorBar',defaultCBar);
defaultCLimits = 'auto';
addParameter(p,'ColorLimits',defaultCLimits);
defaultColorMap = hsv(128);
addParameter(p,'ColorMap',defaultColorMap);
defaultEdgeColor = 'None'; % black
addParameter(p,'EdgeColor',defaultEdgeColor);
defaultFilename = '';
addParameter(p,'Filename',defaultFilename);
defaultHandle = [];
addParameter(p,'Handle',defaultHandle);
defaultMarker = '-';
addParameter(p,'Marker',defaultMarker,@isstr);
defaultPQ = -20:0.125:20;
addParameter(p,'PQ',defaultPQ,@isvector);
defaultStyle = 'advanced';
addParameter(p,'Style',defaultStyle,@isstr);
defaultZLim = [];
addParameter(p,'ZLim',defaultZLim,@isvector);
defaultZString = 'W(q,p)';
addParameter(p,'ZString',defaultZString,@isstr);
parse(p,varargin{:});
c = struct2cell(p.Results);
[cbar,climits,cmap,edgecolor,filename,handle,marker,pq,style,zlimit, ...
    zstring] = c{:};

if strcmp(style,'2D')
    image = true;
else
    image = false;
end

%% Preprocess data and axes
% if narrow
%     p = -6:0.125:6;
%     q = p;
%     nP = length(p);
%     nWF = length(WF);
%     shift = (nWF-nP)/2;
%     WF = WF(shift+1:end-shift,shift+1:end-shift);
% else
%     p = -20:0.125:20;
%     q = p;
% end
[p,q] = deal(pq);
nP = length(p);
nWF = length(WF);
shift = (nWF-nP)/2;
WF = WF(shift+1:end-shift,shift+1:end-shift);
WF = real(WF);

%% Plot data
if ~isempty(handle)
    set(gcf,'currentaxes',handle);
    fig = gcf;
else
    fig = figure;
end

% Prepare figure for export to pdf
% if isempty(handle)
%     formatFigA5(fig);
% end
set(fig,'Color','w');

% Prepare automatic colormap
if strcmp(climits,'auto')
    %cmap = defaultColorMap;
    minVal = min(min(WF));
    maxVal = max(max(WF));
    shift = ceil(abs(minVal)/(maxVal-minVal)*length(cmap));
    cmap = circshift(cmap,shift,1);
end

if strcmp(style,'advanced')
    % Add 3D plot with important details
    colormap(cmap);
    surf(q,p,WF,'EdgeColor',edgecolor,'FaceLighting','gouraud');
    caxis(climits);
    hold on;
    view(3);
    xlim([min(q),max(q)]);
    ylim([min(q),max(q)]);
    if ~isempty(zlimit)
        zlim(zlimit);
    else
        zl = zlim;
        minZl = max(abs(zl));
        zlim([-minZl,zl(2)]);
    end
    if ~isempty(edgecolor)
        camlight('left');
    end
    xlabel('p','FontWeight','bold','FontSize',32);
    ylabel('q','FontWeight','bold','FontSize',32);
    zlabel(zstring,'FontWeight','bold','FontSize',32);
    set(gca,'FontSize',22);
    grid on;

    % Draw image beneath 3D plot and two grid lines
    imgzposition = min(zlim);
    lenq = length(q);
    lenp = length(p);
    surf([min(q) max(q)],[min(p) max(p)],repmat(imgzposition, [2 2]),...
        WF,'facecolor','texture');
    plot3(ones(lenq,1)*0,q,ones(lenq,1)*imgzposition,'Color','black');
    plot3(p,ones(lenp,1)*0,ones(lenp,1)*imgzposition,'Color','black');

    % Draw integral projections along q and p axes
%     searchWF = abs(WF);
%     maxQ = max(searchWF,[],1);
%     [~,indProjP] = max(maxQ);
%     maxP = max(searchWF,[],2);
%     [~,indProjQ] = max(maxP);
    prq = sum(WF,2);
    prq = prq * max(max(WF))/max(prq);
    prp = sum(WF,1);
    prp = prp * max(max(WF))/max(prp);
    plot3(ones(length(p))* max(p),p,prq,marker,'Color','black');
    plot3(q,ones(length(q))* max(p),prp,marker,'Color','black');
    hold off;
    
    % Add colorbar and move it to prevent overlap with xlabel
    if strcmp(cbar,'on')
        colorbar;
    end
else
    colormap(cmap);
    surf(gca,p,q,real(WF),'EdgeColor',edgecolor);
    if ~strcmp(climits,'auto')
        caxis(climits);
    end
    if image
        view(2);
    end
    xlim([min(p),max(p)]);
    ylim([min(q) max(q)]);
    if ~isempty(zlimit)
        zlim(zlimit);
    end
    cb = colorbar;
    if isempty(handle) % adapt figure size only when no figure was prepared
        p = get(gca,'Position');
        set(gca,'OuterPosition',p);
    else
        cb.Label.FontWeight = 'bold';
        cb.Label.FontSize = 24;
        cb.Label.String = zstring;
    end
    title(zstring,'FontSize',32); % old was 32
    xlabel('q','FontWeight','bold','FontSize',32); % old was 32
    ylabel('p','FontWeight','bold','FontSize',32);
    set(gca,'FontSize',22); % old was 22
    set(gca,'DataAspectRatio',[1,1,1]);
end

savefig([filename '.fig']); % does not save empty filename
print([filename '.png'],'-dpng','-r300');

end

