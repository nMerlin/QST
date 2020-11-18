function [H, binsO1, binsO2,r,nTherm, nThermErr,nCoherent,nCohErr,meanN,Coherence] = plotHusimiAndCut(O1,O2,H,binsO1,binsO2,varargin)
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
x = 1/sqrt(2)*sqrt(XAxis.^2+YAxis.^2); % xfit = |alpha|
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
%     % old method: find optimum nTherm and nCoherent, which minimizes the function mini
%     x0 = [nTherm,nCoherent];
%     f = @(x)mini(x, WLength, WRes,H);
%     % set constraints, so nTherm + nCoherent = meanN and both are >= 0
%     A = []; b = []; Aeq = [1,1]; beq = meanN; lb = [0,0]; ub = [Inf,Inf];
%     [x,~] = fmincon(f,x0,A,b,Aeq,beq,lb,ub);
%     nTherm = x(1);
%     nCoherent = x(2); 
%     r = sqrt(2*nCoherent);
%     [ HusFunc,~,~ ] = CreateDisplacedThermalHusimiPhaseAveraged( nTherm, WLength, WRes, r,'PlotOption',false);
%     [~,~,~,~,~,hessian] = fminunc(f,x);
%     err = sqrt(diag(inv(hessian)));
%     nThermErr = err(1);
%     nCohErr = err(2);
%     %The key to the standard errors is the Hessian matrix. 
%     %The variance-covariance-matrix of the coefficients is the inverse
%     % of the Hessian matrix. So the standard errors are the square root of
%     %the values on the diagonal of the inverse Hessian matrix.
%     % The hessian of fminunc is accurate, the one of fmincon is not.
%     % https://de.mathworks.com/matlabcentral/answers/153414-estimator-standard-errors-using-fmincon-portfolio-optimization-context
    
    %% Method with theory function 
    xfit = x(:); %1D
    Hfit = H(:);
    HTheory = fittype('0.5*WRes^2*(pi*(a1+1))^-1 *exp(-(x.^2 + b1)/(a1+1)) .* besseli(0,2*x*sqrt(b1)/(a1+1))','problem','WRes'); 
    %a1 = nTherm; b1 = |alpha0|^2 = nCoherent; 0.5*WRes^2 is for normalization  
    [f,gof,~] = fit(xfit,Hfit,HTheory,'problem',WRes,'StartPoint', [nTherm,nCoherent],'Lower',[0,0] );
    nTherm = f.a1;    
    nCoherent = f.b1;
    level = 2*tcdf(-1,gof.dfe);
    m = confint(f,1-level);    
    nThermErr = m(2,1) - nTherm; 
    nCohErr = m(end,end) - nCoherent; 
end
a = ( nTherm +1)^2 -  nTherm ^2;
Coherence = (1 - exp(-2*nCoherent/a) * besseli(0,2*nCoherent/a ) )/a;
HusFunc = 0.5*WRes^2*(pi*(nTherm+1))^-1 *exp(-(x.^2 + nCoherent)/(nTherm+1)) .* besseli(0,2*x*sqrt(nCoherent)/(nTherm+1));
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

