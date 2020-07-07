function [ HusFunc,xx,yy ] = CreateDisplacedThermalHusimiPhaseAveraged( n, WLength, WRes, r,varargin)
%% Validate and parse input arguments
p = inputParser;
defaultPlotOption = true;
addParameter(p,'PlotOption',defaultPlotOption,@islogical);
parse(p,varargin{:});
c = struct2cell(p.Results);
[plotOption] = c{:};

%%---

MatSize=round(2*(WLength/WRes)+1);
HusFunc=zeros(MatSize,MatSize);

for theta = 0:0.1:2*pi; 
    Q0 = r*cos(theta);
    P0 = r*sin(theta);       
   [ HusFuncPhase,xx,yy ] = CreateDisplacedThermalHusimi( n, Q0,P0,WLength, WRes );
    HusFunc = HusFunc + HusFuncPhase;
end;
HusFunc=HusFunc./sum(sum(HusFunc));

%% Plot stuff
if plotOption
    subplot(3,1,1);
    imagesc(xx(1,:),yy(:,1),HusFunc); axis on; colormap hot; hold on;
    shading('flat');
    xlabel('q');
    ylabel('p');
    title('H(q,p)');
    pbaspect([1 1 1]); % ensure that a circle doesn't look elliptic

    %%
    ax2 = subplot(3,1,2);
    Hcut = HusFunc(round(MatSize/2),:);
    Hcut = Hcut/max(Hcut);
    plot(xx(1,:),Hcut,'DisplayName','measured Husimi');
    xlabel('q');
    ylabel('normalized cut histogram');
    legend(ax2,'location','best');
    ax2.XTick = (min(ax2.XLim):1:max(ax2.XLim));

    %%
    ax3 = subplot(3,1,3);
    Hcut = HusFunc(:,round(MatSize/2));
    Hcut = Hcut/max(Hcut);

    plot(yy(:,1),Hcut,'DisplayName','measured Husimi');
    xlabel('p');
    ylabel('normalized cut histogram');
    legend(ax3,'location','best');
    ax3.XTick = (min(ax3.XLim):1:max(ax3.XLim));
    % savefig([filename '-nbins-' num2str(nBinsA) '-Husimi.fig']);
% print([filename '-nbins-' num2str(nBinsA) '-Husimi.png'],'-dpng');
end






end
