function [H, binsO1, binsO2,r,nTherm, nThermErr,nCoherent,nCohErr,meanN,Coherence,CoherenceErr,poissonErrors] = plotHusimiAndCut(O1,O2,H,binsO1,binsO2,varargin)
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
defaultFitMethod = 'NLSQ-LAR'; %this is the most stable method 
addParameter(p,'FitMethod',defaultFitMethod,@isstr);
defaultiSelect = []; % Indices of selected points shown in green. 
addParameter(p,'iSelect',defaultiSelect,@isvector);
defaultShowLegend = true;
addParameter(p,'ShowLegend',defaultShowLegend,@islogical);
defaultPlotOption = true;
addParameter(p,'PlotOption',defaultPlotOption,@islogical);
defaultMonteCarloError = false; %In this case, this function is only used to 
%get errors from MonteCarlo error estimation in HusimiSeriesFromQuadratures.m
addParameter(p,'MonteCarloError',defaultMonteCarloError,@islogical);
defaultFitFunction = []; %the fitFunction used to fit. If it is empty, a new fit function will be created. 
addParameter(p,'FitFunction',defaultFitFunction);
parse(p,varargin{:});
c = struct2cell(p.Results);
[filename,fitFunction,fitMethod,iSelect,limits,monteCarloError,nBinsA,nBinsB,plotOption,showLegend] = c{:};

%% compute Husimi, if it is not given
%(not needed when we already have Husimi for Monte Carlo Error estimation)
if ~monteCarloError
    O1 = O1*sqrt(2);
    O2 = O2*sqrt(2);
%     For the Husimi function, the quadratures are scaled so they correspond to
    % those quadratures before the Husimi beam splitter
    if H == 0
        [H, binsO1, binsO2] = histogram2D(O1,O2,'nBinsA',nBinsA,'nBinsB',nBinsB);
    end;
    H=H./(sum(sum(H)));
    totalNumber = length(O1);
    poissonErrors = sqrt(H.*(1-H)/totalNumber);
    % a matrix with the counting error on each phase space bin. 
    poissonErrorsCut = poissonErrors(round(length(binsO2)/2),:);
else
    poissonErrors = 0;
    poissonErrorsCut = 0;
end

%% get mean total photon number from averaging over Husimi function 
[XAxis,YAxis]=meshgrid(binsO1,binsO2);
x = 1/sqrt(2)*sqrt(XAxis.^2+YAxis.^2); % = |alpha|
FullQx2=sum(sum(XAxis.^2.*H));
FullQy2=sum(sum(YAxis.^2.*H));
meanN=0.5*(FullQx2+FullQy2) -1;
%% get cut of H along q axis
Hcut = H(round(length(binsO2)/2),:);

%% Get Theory
WLength = max(binsO1);
WRes = (2*WLength)/(nBinsA - 1);
%get starting values for r and nTherm 
ys = transpose(csaps(binsO1,Hcut,0.6,binsO1));
[~,I] = max(ys);
r = abs(binsO1(I));
nCoherent = 0.5*r^2;
nTherm = meanN-nCoherent; % thermal state photon number
xfit = x(:); %1D
Hfit = H(:);
    
if strcmp(fitMethod,'NLSQ-LAR')   
    if isempty(fitFunction)
        fitFunction = fittype('0.5*WRes^2*(pi*(a1+1))^-1 *exp(-(x.^2 + b1)/(a1+1)) .* besseli(0,2*x*sqrt(b1)/(a1+1))','problem','WRes'); 
    end
    %a1 = nTherm; b1 = |alpha0|^2 = nCoherent; 0.5*WRes^2 is for normalization  
    [f,gof,~] = fit(xfit,Hfit,fitFunction,'problem',WRes,'StartPoint', [nTherm,nCoherent],'Lower',[0,0],'Robust','LAR' );
    nTherm = f.a1;    
    nCoherent = f.b1;
    if ~monteCarloError
        [se]= getStandardErrorsFromFit(f,gof,'method1');   
        nThermErr = se(1);
        nCohErr = se(2);
    else
         [nThermErr, nCohErr] = deal(0);
    end
elseif strcmp(fitMethod,'NLSQ-Bisquare')
    if isempty(fitFunction)
        fitFunction = fittype('0.5*WRes^2*(pi*(a1+1))^-1 *exp(-(x.^2 + b1)/(a1+1)) .* besseli(0,2*x*sqrt(b1)/(a1+1))','problem','WRes'); 
    end
    %a1 = nTherm; b1 = |alpha0|^2 = nCoherent; 0.5*WRes^2 is for normalization      
    [f,gof,~] = fit(xfit,Hfit,fitFunction,'problem',WRes,'StartPoint', [nTherm,nCoherent],'Lower',[0,0],'Robust','Bisquare' );
    nTherm = f.a1;    
    nCoherent = f.b1;
    [se]= getStandardErrorsFromFit(f,gof,'method1');   
    nThermErr = se(1);
    nCohErr = se(2);
elseif strcmp(fitMethod,'NLSQ-Off')
    if isempty(fitFunction)
        fitFunction = fittype('0.5*WRes^2*(pi*(a1+1))^-1 *exp(-(x.^2 + b1)/(a1+1)) .* besseli(0,2*x*sqrt(b1)/(a1+1))','problem','WRes'); 
    end 
    %a1 = nTherm; b1 = |alpha0|^2 = nCoherent; 0.5*WRes^2 is for normalization      
    [f,gof,~] = fit(xfit,Hfit,fitFunction,'problem',WRes,'StartPoint', [nTherm,nCoherent],'Lower',[0,0],'Robust','Off' );
    nTherm = f.a1;    
    nCoherent = f.b1;
    [se]= getStandardErrorsFromFit(f,gof,'method1');   
    nThermErr = se(1);
    nCohErr = se(2);
elseif strcmp(fitMethod,'fmincon')
    % old method: find optimum nTherm and nCoherent, which minimizes the function f
    p0 = [nTherm,nCoherent];  
    %f = @(p)mini(p, WLength, WRes,H); 
    f = @(p) sum(abs(Hfit - 0.5*WRes^2*(pi*(p(1)+1))^-1 *exp(-(xfit.^2 + p(2))/(p(1)+1)) .* besseli(0,2*xfit*sqrt(p(2))/(p(1)+1))));
    % set constraints, so nTherm + nCoherent = meanN and both are >= 0
    A = []; b = []; Aeq = [1,1]; beq = meanN; lb = [0,0]; ub = [Inf,Inf];
    [p,~] = fmincon(f,p0,A,b,Aeq,beq,lb,ub);
    nTherm = p(1);
    nCoherent = p(2); 
    [~,~,~,~,~,hessian] = fminunc(f,p);
    err = sqrt(diag(inv(hessian)));
    nThermErr = err(1);
    nCohErr = err(2);
%     %The key to the standard errors is the Hessian matrix. 
%     %The variance-covariance-matrix of the coefficients is the inverse
%     % of the Hessian matrix. So the standard errors are the square root of
%     %the values on the diagonal of the inverse Hessian matrix.
%     % The hessian of fminunc is accurate, the one of fmincon is not.
%     % https://de.mathworks.com/matlabcentral/answers/153414-estimator-standard-errors-using-fmincon-portfolio-optimization-context  
else    
    [nThermErr, nCohErr] = deal(0);
end

Coherence = coherencePDTS(nTherm,nCoherent);
if monteCarloError
    CoherenceErr = 0;
else
    [~, CoherenceErr,~, ~] = error_propagation( @(nTherm,nCoherent) coherencePDTS(nTherm,nCoherent), ...
        nTherm, nCoherent, nThermErr, nCohErr );
    CoherenceErr(isnan(CoherenceErr)) = 0;
    HusFunc = 0.5*WRes^2*(pi*(nTherm+1))^-1 *exp(-(x.^2 + nCoherent)/(nTherm+1)) .* besseli(0,2*x*sqrt(nCoherent)/(nTherm+1));
    theoryHFCut = HusFunc(round(nBinsA/2),:);
end

%% plot
if plotOption
    %% plot Husimi function
    pcolor(binsO1,binsO2,H); shading 'flat'; axis on; colormap hot; colorbar; hold on;
    scatter(O1(iSelect),O2(iSelect),'.g','DisplayName','Postselection'); 
    % plot(binsO1,Hcut/max(Hcut)*0.5*max(binsO2)-max(binsO2),'w','Linewidth',2,'DisplayName','Cut');
    % hold on;
    % plot(binsO1,theoryHFCut/max(theoryHFCut)*0.5*max(binsO2)-max(binsO2),'r','Linewidth',2,'DisplayName',...
    %     ['Theory, n_{Th} = ' num2str(nTherm,'%.1f') ', n_{Coh} = ' num2str(nCoherent,'%.1f') ]);
    %set(gca,'XLim',limits,'YLim',limits);
    xlabel('q');
    ylabel('p');
    pbaspect([1 1 1]);
    graphicsSettings;grid;
    ax = gca;
    set(ax,'FontSize',36,'FontName','Arial', 'TickDir','out');
    if showLegend
        legend('location','bestoutside','Fontsize',10);
    end
    savefig([filename '-nbins-' num2str(nBinsA) '-fitMethod-' fitMethod '-Husimi.fig']);
    print([filename '-nbins-' num2str(nBinsA) '-fitMethod-' fitMethod '-Husimi.png'],'-dpng','-r700');
    clf();


    %%  plot Cut through Husimi function 
    line = shadedErrorBar(binsO1,Hcut,poissonErrorsCut);
    line.DisplayName= 'Data';
    hold on;
    plot(binsO1,theoryHFCut,'r','Linewidth',2,'DisplayName',...
        ['Theory, n_{Th} = ' num2str(nTherm,'%.0f') ', n_{Coh} = ' num2str(nCoherent,'%.0f') ]);
    xlabel('q');
    ylabel('Q(q,p)');
    %legend('location','southwest');
    graphicsSettings;
    ax = gca;
    set(ax,'FontSize',36,'FontName','Arial');
    savefig([filename '-nbins-' num2str(nBinsA) '-fitMethod-' fitMethod '-Cut.fig']);
    print([filename '-nbins-' num2str(nBinsA) '-fitMethod-' fitMethod '-Cut.png'],'-dpng','-r700');
end



end

