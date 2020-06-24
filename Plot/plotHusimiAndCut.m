function plotHusimiAndCut(O1,O2,X1,X2,varargin)
%Plot histogram of Husimi-Q function and cut along q axis
%% Validate and parse input arguments
p = inputParser;
defaultFilename = '';
addParameter(p,'Filename',defaultFilename);
defaultLimits = [-7,7];
addParameter(p,'Limits',defaultLimits);
defaultnBinsA = 1000;
addParameter(p,'nBinsA',defaultnBinsA,@isnumeric);
defaultnBinsB = 1000;
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
Hcut = Hcut/sum(Hcut);
[n1,n2,~] = nPhotons(X1,X2,X2);
meanN = mean([n1 n2]);

cohHF = cohHusimi(binsO1,binsO2,meanN);
cohHFMirror = flip(cohHF,2);
cohHFCut = cohHF(round(length(binsO2)/2),:) + cohHFMirror(round(length(binsO2)/2),:);
cohHFCut = cohHFCut/sum(cohHFCut);

plot(binsO1,Hcut,'DisplayName','measured Husimi');
hold on;
plot(binsO1,cohHFCut,'Linewidth',2,'DisplayName','theory coherent Husimi');
xlabel('q');
ylabel('normalized cut histogram');
legend(ax2,'location','best');
ax2.XTick = (min(ax2.XLim):1:max(ax2.XLim));

%%
ax3 = subplot(3,1,3);
Hcut = H(:,round(length(binsO1)/2));
Hcut = Hcut/sum(Hcut);

plot(binsO2,Hcut,'DisplayName','measured Husimi');
hold on;
plot(binsO2,cohHFCut,'Linewidth',2,'DisplayName','theory coherent Husimi');
xlabel('p');
ylabel('normalized cut histogram');
legend(ax3,'location','best');
ax3.XTick = (min(ax3.XLim):1:max(ax3.XLim));

savefig([filename '-nbins-' num2str(nBinsA) '-Husimi.fig']);
print([filename '-nbins-' num2str(nBinsA) '-Husimi.png'],'-dpng');

end

