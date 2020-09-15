function [H, binsO1, binsO2,r,nTherm,nCoherent,meanN] = plotHusimiAndCut(O1,O2,H,binsO1,binsO2,varargin)
%Plot histogram of Husimi-Q function and cut along q axis
% Either make it new from the orthogonal postselection quadratures O1,O2.
% Then set H, binsO1,binsO2 = 0.
% Or make it from already obtained H, binsO1,binsO2.
% The photon numbers correspond to before the Husimi beam splitter.
%% Validate and parse input arguments
p = inputParser;
defaultFilename = '';
addParameter(p,'Filename',defaultFilename);
defaultLimits = [-20,20];
addParameter(p,'Limits',defaultLimits);
defaultnBinsA = 100;
addParameter(p,'nBinsA',defaultnBinsA,@isnumeric);
defaultnBinsB = 100;
addParameter(p,'nBinsB',defaultnBinsB,@isnumeric);
defaultMakeFit = true;
addParameter(p,'MakeFit',defaultMakeFit,@islogical);
parse(p,varargin{:});
c = struct2cell(p.Results);
[filename,limits,makeFit,nBinsA,nBinsB] = c{:};

%% For the Husimi function, the quadratures are scaled so they correspond to
% those quadratures before the Husimi beam splitter
O1 = O1*sqrt(2);
O2 = O2*sqrt(2);
%% plot Husimi function
subplot(3,1,1);

if H == 0
    plotHusimi(O1,O2,[],'Limits',limits,'nBinsA',nBinsA,'nBinsB',nBinsB); % husimi=[O1,O2,iSelect]
    %%
    [H, binsO1, binsO2] = histogram2D(O1,O2,'nBinsA',nBinsA,'nBinsB',nBinsB);
    H=H./(sum(sum(H)));
else
    H=H./(sum(sum(H)));
    imagesc(binsO1,binsO2,H); axis on; colormap hot; hold on;
    set(gca,'XLim',limits,'YLim',limits);
    xlabel('q');
    ylabel('p');
    title('H(q,p)');
    pbaspect([1 1 1]); 
end

%% get mean total photon number from averaging over Husimi function 
[XAxis,YAxis]=meshgrid(binsO1,binsO2);
FullQx2=sum(sum(XAxis.^2.*H));
FullQy2=sum(sum(YAxis.^2.*H));
meanN=0.5*(FullQx2+FullQy2) -1;
%% get cut of H along q axis
Hcut = H(round(length(binsO2)/2),:);
Hcut = Hcut/max(Hcut);


%% Get Theory
WLength = max(binsO1);
WRes = (2*WLength)/(nBinsA - 1);
%get starting values for r and nTherm 
ys = transpose(csaps(binsO1,Hcut,0.6,binsO1));
[~,I] = max(ys);
r = abs(binsO1(I));
nCoherent = 0.5*r^2;
nTherm = meanN-nCoherent; % thermal state photon number

if makeFit
    % find optimum nTherm and nCoherent, which minimizes the function mini
    x0 = [nTherm,nCoherent];
    f = @(x)mini(x, WLength, WRes,H);
    % set constraints, so nTherm + nCoherent = meanN and both are >= 0
    A = []; b = []; Aeq = [1,1]; beq = meanN; lb = [0,0]; ub = [Inf,Inf];
    [x,~] = fmincon(f,x0,A,b,Aeq,beq,lb,ub);
    nTherm = x(1);
    nCoherent = x(2); 
    r = sqrt(2*nCoherent);
end

[ HusFunc,~,~ ] = CreateDisplacedThermalHusimiPhaseAveraged( nTherm, WLength, WRes, r,'PlotOption',false);
theoryHFCut = HusFunc(round(nBinsA/2),:);
theoryHFCut = theoryHFCut/max(theoryHFCut);

%% plot stuff 
ax2 = subplot(3,1,2);
plot(binsO1,Hcut,'DisplayName','measured Husimi');
hold on;
plot(binsO1,theoryHFCut,'Linewidth',2,'DisplayName',...
    ['theory displaced thermal Husimi, r = ' num2str(r) ', n_{Th} = ' num2str(nTherm) ', n_{Coh} = ' num2str(nCoherent) ]);
xlabel('q');
ylabel('normalized cut histogram');
legend(ax2,'location','southwest');
ax2.XTick = (min(ax2.XLim):1:max(ax2.XLim));

%%
ax3 = subplot(3,1,3);
Hcut = H(:,round(length(binsO1)/2));
Hcut = Hcut/max(Hcut);

plot(binsO2,Hcut,'DisplayName','measured Husimi');
hold on;
plot(binsO2,theoryHFCut,'Linewidth',2,'DisplayName',...
    ['theory displaced thermal Husimi, r = ' num2str(r) ', n_{Th} = ' num2str(nTherm) ', n_{Coh} = ' num2str(nCoherent) ]);
xlabel('p');
ylabel('normalized cut histogram');
legend(ax3,'location','southwest');
ax3.XTick = (min(ax3.XLim):1:max(ax3.XLim));

savefig([filename '-nbins-' num2str(nBinsA) '-Husimi.fig']);
print([filename '-nbins-' num2str(nBinsA) '-Husimi.png'],'-dpng');

end

