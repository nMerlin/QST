function plotHusimiAndCut(O1,O2,X1,X2,r,varargin)
%Plot histogram of Husimi-Q function and cut along q axis
% O1, O2: orthogonal postselection quadratures 
% r: radius of coherent state
%% Validate and parse input arguments
p = inputParser;
defaultFilename = '';
addParameter(p,'Filename',defaultFilename);
defaultLimits = [-10,10];
addParameter(p,'Limits',defaultLimits);
defaultnBinsA = 100;
addParameter(p,'nBinsA',defaultnBinsA,@isnumeric);
defaultnBinsB = 100;
addParameter(p,'nBinsB',defaultnBinsB,@isnumeric);
parse(p,varargin{:});
c = struct2cell(p.Results);
[filename,limits,nBinsA,nBinsB] = c{:};



%%
subplot(3,1,1);
plotHusimi(O1,O2,[],'Limits',limits,'nBinsA',nBinsA,'nBinsB',nBinsB); % husimi=[O1,O2,iSelect]
 
%%
ax2 = subplot(3,1,2);
[H, binsO1, binsO2] = histogram2D(O1,O2,'nBinsA',nBinsA,'nBinsB',nBinsB);
Hcut = H(round(length(binsO2)/2),:);
Hcut = Hcut/max(Hcut);
[n1,n2,~] = nPhotons(X1,X2,X2);
meanN = mean([n1 n2]);

%% Get Theory
WLength = max(binsO1);
WRes = (2*WLength)/(nBinsA - 1);
nTherm = meanN-0.5*r^2; % thermal state photon number
[ HusFunc,~,~ ] = CreateDisplacedThermalHusimiPhaseAveraged( nTherm, WLength, WRes, r,'PlotOption',false);
theoryHFCut = HusFunc(round(nBinsA/2),:);
theoryHFCut = theoryHFCut/max(theoryHFCut);

% cohHF = cohHusimi(binsO1,binsO2,meanN);
% cohHFMirror = flip(cohHF,2);
% theoryHFCut = cohHF(round(length(binsO2)/2),:) + cohHFMirror(round(length(binsO2)/2),:);
% theoryHFCut = theoryHFCut/max(theoryHFCut);

plot(binsO1,Hcut,'DisplayName','measured Husimi');
hold on;
plot(binsO1,theoryHFCut,'Linewidth',2,'DisplayName',...
    ['theory displaced thermal Husimi, r = ' num2str(r) ', n_{Th} = ' num2str(nTherm)]);
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
    ['theory displaced thermal Husimi, r = ' num2str(r) ', n_{Th} = ' num2str(nTherm)]);
xlabel('p');
ylabel('normalized cut histogram');
legend(ax3,'location','southwest');
ax3.XTick = (min(ax3.XLim):1:max(ax3.XLim));

savefig([filename '-nbins-' num2str(nBinsA) '-Husimi.fig']);
print([filename '-nbins-' num2str(nBinsA) '-Husimi.png'],'-dpng');

end

