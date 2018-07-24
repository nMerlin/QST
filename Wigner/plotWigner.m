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
defaultFilename = '';
addParameter(p,'Filename',defaultFilename);
defaultHandle = [];
addParameter(p,'Handle',defaultHandle);
defaultPQ = -20:0.125:20;
addParameter(p,'PQ',defaultPQ,@isvector);
defaultStyle = '3D';
addParameter(p,'Style',defaultStyle,@isstr);
defaultZLim = [];
addParameter(p,'ZLim',defaultZLim,@isvector);
parse(p,varargin{:});
c = struct2cell(p.Results);
[filename,handle,pq,style,zlimit] = c{:};

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

%% Plot data
if ~isempty(handle)
    set(gcf,'currentaxes',handle);
    fig = gcf;
else
    fig = figure;
end

if strcmp(style,'advanced')
    % Prepare figure for export to pdf
    hold on;
    formatFigA5(fig);
    set(fig,'Color','w');

    % Add 3D plot with important details
    surf(q,p,WF);
    if ~isempty(zlimit)
        zlim(zlimit);
    else
        zl = zlim;
        minZl = max(abs(zl));
        zlim([-minZl,zl(2)]);
    end
    xlabel('q','FontWeight','bold','FontSize',32);
    ylabel('p','FontWeight','bold','FontSize',32);
    zlabel('W(q,p)','FontWeight','bold','FontSize',32);
    set(gca,'FontSize',22);
    grid on;
    view(3);

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
    plot3(ones(length(p))*5,p,prq,'Color','black');
    plot3(q,ones(length(q))*5,prp,'Color','black');
    hold off;
else
    h = surf(gca,p,q,real(WF));
    set(h,'LineStyle','none');
    if image
        view(2);
    end
    xlim([min(p),max(p)]);
    ylim([min(q) max(q)]);
    zlim(zlimit);
end

if ~isempty(filename)
    [path,name] = fileparts(filename);
    filenamePdf = [path,name,'.pdf'];
    export_fig(sprintf('%s',filenamePdf),'-pdf','-painters');
end

end

