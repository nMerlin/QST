function plotPFunctionResults(listOfParams,varargin)
%PLOTSERIESPOSTSELECTIONS

%% Validate and parse input arguments
p = inputParser;
defaultXUnit = 'ps';
addParameter(p,'XUnit',defaultXUnit,@isstr);
defaultFitType = 'gauss';
addParameter(p,'FitType',defaultFitType,@isstr);
defaultVaryAPS = false;
addParameter(p,'VaryAPS',defaultVaryAPS,@islogical);
defaultRemoveModulation = false;
addParameter(p,'RemoveModulation',defaultRemoveModulation,@islogical);
defaultRange = 0.3;
addParameter(p,'Range',defaultRange,@isvector);
defaultZeroDelay = 0;
addParameter(p,'ZeroDelay',defaultZeroDelay,@isnumeric); % in mm
defaultPlotrelative = false; %plots the data relative to the target photon number
addParameter(p,'Plotrelative',defaultPlotrelative,@islogical);
defaultLogplot = false; %makes yscale of plot logarithmic
addParameter(p,'Logplot',defaultLogplot,@islogical);
defaultMaxQuad = 20; %max abs value of the new quadrature coordinates 
addParameter(p,'MaxQuad',defaultMaxQuad,@isnumeric);
defaultMaxX= 20; %max abs value of the old quadrature coordinates 
addParameter(p,'MaxX',defaultMaxX,@isnumeric);
defaultPhiStep= 0.1; %step size of the phase grid
addParameter(p,'PhiStep',defaultPhiStep,@isnumeric);
defaultRvalue= 0.7; %filter Parameter R
addParameter(p,'Rvalue',defaultRvalue,@isnumeric);
defaultXStep= 1; %step size of the old quadrature coordinates grid
addParameter(p,'XStep',defaultXStep,@isnumeric);
defaultResolution= 1; %step size of the new quadrature coordinates grid
addParameter(p,'Resolution',defaultResolution,@isnumeric);
defaultNorm= 1; %normalization of the vacuum standard deviation of quadratures. Is 1 for P functions
addParameter(p,'Norm',defaultNorm,@isnumeric);
parse(p,varargin{:});
c = struct2cell(p.Results);
[fitType,logplot,maxQuad,maxX,norm,phiStep,plotrelative,range,remMod,res,rvalue,varyAPS,XStep,xUnit,zeroDelay] = c{:};

%% Create folder 'figures-fig'
folder = ['Pfunction-figures-results-R-' num2str(rvalue) ...
        '-maxQuad-' num2str(maxQuad) '-Resolution-' num2str(res) '-maxX-' num2str(maxX) '-Xstep-' num2str(XStep) '-phiStep-' num2str(phiStep)];
if ~exist([pwd folder],'dir')
    mkdir(folder);
end

% Constants
figurepath = [folder '/'];
load('Photonnumbers.mat','nPsFast','nPsSlow','nTg'); %loads the photon numbers of each channel without postselection,...
%which was created with plotSeriesPostselections
%% Gather data
[delay,Yr,Yt,PhotonNrs,meanPh,meanR,varPh,varR,Rerr,phErr,varRerr,varPhErr,PhotonNrErrs] = deal([]);
for iParams = 1:length(listOfParams)
    selParams = listOfParams(iParams);
    selStr = selParamsToStr(selParams);
    foldername = ['Pfunctionplots-',selStr,'-remMod-',...
            num2str(remMod),'-range-',num2str(range),'-varyAPS-',num2str(varyAPS),'-R-' num2str(rvalue) ...
        '-maxQuad-' num2str(maxQuad) '-Resolution-' num2str(res) '-maxX-' num2str(maxX) '-Xstep-' num2str(XStep) '-phiStep-' num2str(phiStep)];
    load([foldername '\Pfunctionresults.mat'],'Delay','DelayMm','Pmax','sigmaPmax','meanPhase','meanPhaseBinned','varPhase',...
    'varPhaseBinned','g1','sigNeg','meanAmp','meanAmpBinned','varAmp','varAmpBinned','PhotonNr','PhotonNrBinned',...
    'phaseErr','ampErr','varPhaseErr','varAmpErr','PhotonNrErr');
    delay(iParams,:) = Delay; 
    H = length(delay);
    Radii = ones(H,1) * selParams.Position(1);
    Thicknesses = ones(H,1) * selParams.Position(2);
    Yr(iParams,:) = Radii;
    Yt(iParams,:) = Thicknesses;
    PhotonNrs(iParams,:) = PhotonNr;
    meanPh(iParams,:) = meanPhase;
    meanR(iParams,:) = meanAmp; 
    varPh(iParams,:) = varPhase;
    varR(iParams,:) = varAmp;
    Rerr(iParams,:) = ampErr;
    phErr(iParams,:) = phaseErr;
    varRerr(iParams,:) = varAmpErr;
    varPhErr(iParams,:) = varPhaseErr;
    PhotonNrErrs(iParams,:) = PhotonNrErr;   
end
[~,I] = sort(delay(:,1)); % Sort for Radii
delay = delay(I,:); %delays
Yr = Yr(I,:); %Radii
Yt = Yt(I,:); %thicknesses
if norm == 1
    Yr = Yr*sqrt(2);
    Yt = Yt*sqrt(2);
end  
PhotonNrs = PhotonNrs(I,:);
meanPh = meanPh(I,:);
meanR = meanR(I,:);varPh = varPh(I,:);varR = varR(I,:);
Rerr = Rerr(I,:); phErr = phErr(I,:); varRerr = varRerr(I,:); varPhErr = varPhErr(I,:); PhotonNrErrs = PhotonNrErrs(I,:);
[fitTau,tauErr,fitPeak,fitPeakErr,sse,rsquare,adjrsquare,rmse,pa1,pa2,pa3,pa4,pa1Err,pa2Err,pa3Err,pa4Err] = deal(zeros(length(I),1));
%fitTau: coherencetime; pa1 etc: fit parameters 
if plotrelative
    meanR = meanR./sqrt(nTg);
    PhotonNrs = PhotonNrs./nTg;
end

if varyAPS
    sel = 'k'; %the selection radius or variable
else
    sel = 'r';
end

%% Plot
typestrVector = {'R','discN','varR','meanPh','varPh','RLog'};
for typeI = 1:length(typestrVector)
    plotStuff(cell2mat(typestrVector(typeI)))
end

    function [] = plotStuff(typestr)
%% Create figure
fig = figure;
filename = [figurepath,typestr,'-',fitType,'-remMod-',...
            num2str(remMod),'-range-',num2str(min(range)),...
            '-',num2str(max(range)),'-varyAPS-',num2str(varyAPS),...
            '-plotrelative-' num2str(plotrelative),...
            '-logplot-',num2str(logplot),'.fig'];
%formatFigA5(fig);
figure(1);
ax = gca;
for i = 1:length(I)
    switch typestr
        case 'R'
            ys = meanR(i,:);
            yErr = Rerr(i,:);
            if plotrelative
                ylabel(ax,'$mean amplitude <r>/\sqrt{n_{Tg}}$','Interpreter','latex'); 
            else            
                ylabel(ax,'Mean amplitude <r>'); 
            end
        case 'RLog'
            ax.YScale = 'log';
            ys = meanR(i,:);
            yErr = Rerr(i,:);
            if plotrelative
                ylabel(ax,'$mean amplitude <r>/\sqrt{n_{Tg}}$','Interpreter','latex'); 
            else            
                ylabel(ax,'mean amplitude <r>'); 
            end    
        case 'meanPh'
            ys = meanPh(i,:);
            yErr = phErr(i,:);
            ylabel(ax,'<\phi> (rad)'); 
        case 'varPh'
            ys = varPh(i,:); 
            yErr = varPhErr(i,:);
            ylabel(ax,'Variance(\phi) (rad)');
        case 'varR'    
            ys = varR(i,:);
            yErr = varRerr(i,:);
            ylabel('Variance(r)');      
        case 'discN'
            ys = PhotonNrs(i,:);
            yErr = PhotonNrErrs(i,:);
            if plotrelative
                ylabel(ax,'$n/n_{Tg}$','Interpreter','latex'); 
            else            
                ylabel(ax,'n'); 
            end
    end %typestr     
    
    shadedErrorBar(delay(i,:),ys,yErr,'lineProps',{'o-','DisplayName',[sel ' = ' num2str(Yr(i,1),2) ', t = ' num2str(Yt(i,1),2)]});
    %plot(ax,delay(i,:),ys,'o-','DisplayName',[sel ' = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
    hold on;
    x = delay(i,:);
    if ~any(~isnan(ys)) %if there are only nans
        ys = zeros(1,H);
    end
%             ys = transpose(csaps(x,ys,0.00001,x));  
%             ys = ys';
    switch fitType
        case 'gauss'
            g = fittype(@(a,b,c,d,x) a*exp(-pi/2*((x-b)/c).^2)+d);                   
            b0 = zeroDelay;
            c0 = 0.5*max(x);
            d0 = mean(ys(end-5:end)); % saturation value
            a0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay))) - d0;
            [result,gof,~] = fit(x',ys',g,'StartPoint',[a0 b0 c0 d0]); 
            [se]= getStandardErrorsFromFit(result,gof,'method1'); 
            pa1(i) = result.a;
            pa2(i) = result.b;
            pa3(i) = result.c;
            pa4(i) = result.d;
            pa1Err(i) = se(1);
            pa2Err(i) = se(2);
            pa3Err(i) = se(3);
            pa4Err(i) = se(4);
            fitTau(i) = result.c;            
            tauErr(i) = se(3);
            fitPeak(i) = result.a + result.d; 
            fitPeakErr(i) = sqrt(se(1)^2 + se(4)^2);           
            plot(ax,min(delay(i,:)):1:max(delay(i,:)),result(min(delay(i,:)):1:max(delay(i,:))),'r','DisplayName','');            
         case 'gauss2'
            g = fittype(@(a,c,x) a*exp(-pi/2*((x)/c).^2));                   
            c0 = 0.5*max(x);
            a0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay)));
            [result,gof,~] = fit(x',ys',g,'StartPoint',[a0 c0]); 
            [se]= getStandardErrorsFromFit(result,gof,'method1'); 
            pa1(i) = result.a;          
            pa2(i) = result.c;
            pa1Err(i) = se(1);
            pa2Err(i) = se(2);
            fitTau(i) = result.c;            
            tauErr(i) = se(2);
            fitPeak(i) = result.a; 
            fitPeakErr(i) = se(1);           
            plot(ax,min(delay(i,:)):1:max(delay(i,:)),result(min(delay(i,:)):1:max(delay(i,:))),'r','DisplayName','');
        case 'exponential'
            ex = fittype(@(m,c,x) m - x/c);                   
            b0 = zeroDelay;
            c0 = 0.5*max(x);                
            d= mean(ys(end-5:end));
            ys = ys-d; 
            a0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay)));
            sig = sign(a0);
            ys = abs(ys);
            a0 = abs(a0);
            m0 = log(a0) + b0/c0;
            [result,gof,~] = fit(abs(x)',log(ys)',ex,'StartPoint',[m0 c0]);             
            fitPeak(i) = sig*exp(result(b0))+d;              
            [se]= getStandardErrorsFromFit(result,gof,'method1'); 
            pa1(i) = result.m;
            pa2(i) = result.c;
            pa1Err(i) = se(1);
            pa2Err(i) = se(2);
            fitTau(i) = result.c;            
            tauErr(i) = se(2);
            fitPeak(i) = sig*exp(result(b0))+d;            
            plot(ax,min(delay(i,:)):1:max(delay(i,:)),sig*exp(result(abs(min(delay(i,:)):1:max(delay(i,:)))))+d,'r','DisplayName',''); 
        case 'exp2'
            g = fittype(@(a,b,c,d,x) a*exp(-((x-b)/c)) +d );                   
            b0 = zeroDelay;
            c0 = 0.5*max(x);
            d0 = mean(ys(end-5:end)); % saturation value
            a0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay))) - d0;
            [result,gof,~] = fit(abs(x)',ys',g,'StartPoint',[a0 b0 c0 d0]); 
            [se]= getStandardErrorsFromFit(result,gof,'method1'); 
            pa1(i) = result.a;
            pa2(i) = result.b;
            pa3(i) = result.c;
            pa4(i) = result.d;
            pa1Err(i) = se(1);
            pa2Err(i) = se(2);
            pa3Err(i) = se(3);
            pa4Err(i) = se(4);
            fitTau(i) = result.c;            
            tauErr(i) = se(3);            
            fitPeak(i) = result.a + result.d; 
            fitPeakErr(i) = sqrt(se(1)^2 + se(4)^2); 
            plot(ax,min(delay(i,:)):1:max(delay(i,:)),result(abs(min(delay(i,:)):1:max(delay(i,:)))),'r','DisplayName','');
        case 'exp3'
            g = fittype(@(a,b,c,x) a*exp(-((x-b)/c))  );                   
            b0 = zeroDelay;
            c0 = 0.5*max(x);
            a0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay)));
            [result,gof,~] = fit(abs(x)',ys',g,'StartPoint',[a0 b0 c0]); 
            [se]= getStandardErrorsFromFit(result,gof,'method1'); 
            pa1(i) = result.a;
            pa2(i) = result.b;
            pa3(i) = result.c;
            pa1Err(i) = se(1);
            pa2Err(i) = se(2);
            pa3Err(i) = se(3);
            fitTau(i) = result.c;            
            tauErr(i) = se(3);            
            fitPeak(i) = result.a; 
            fitPeakErr(i) = se(1);  
            plot(ax,min(delay(i,:)):1:max(delay(i,:)),result(abs(min(delay(i,:)):1:max(delay(i,:)))),'r','DisplayName','');
        case 'sech'
            se = fittype(@(a,b,c,d,x) a./(cosh(-pi/2*(x-b)/c)).^2+d);                   
            b0 = zeroDelay;
            c0 = 0.5*max(x);
            d0 = mean(ys(end-5:end)); % saturation value
            a0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay))) - d0;
            [result,gof,~] = fit(x',ys',se,'StartPoint',[a0 b0 c0 d0]); 
            [se]= getStandardErrorsFromFit(result,gof,'method1'); 
            pa1(i) = result.a;
            pa2(i) = result.b;
            pa3(i) = result.c;
            pa4(i) = result.d;
            pa1Err(i) = se(1);
            pa2Err(i) = se(2);
            pa3Err(i) = se(3);
            pa4Err(i) = se(4);
            fitTau(i) = result.c;            
            tauErr(i) = se(3);            
            fitPeak(i) = result.a + result.d; 
            fitPeakErr(i) = sqrt(se(1)^2 + se(4)^2);  
            plot(ax,min(delay(i,:)):1:max(delay(i,:)),result(min(delay(i,:)):1:max(delay(i,:))),'r','DisplayName','');
        case 'power-law'
            g = fittype(@(a,alpha,x) a*abs(x).^(-alpha));                   
            b0 = zeroDelay;
            alpha0 = 0.5;
            a0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay)));
            [result,gof,~] = fit( x(abs(delay(i,:))>= (100+zeroDelay))',...
                 ys(abs(delay(i,:))>= (100+zeroDelay))',g,'StartPoint',[a0 alpha0]); 
                        [se]= getStandardErrorsFromFit(result,gof,'method1'); 
            pa1(i) = result.a;
            pa2(i) = result.alpha;
            pa1Err(i) = se(1);
            pa2Err(i) = se(2);           
            fitPeak(i) = result.a; 
            fitPeakErr(i) = se(1);  
            % For the power law, there is no coherence time, because an
            % integral over it doesnt converge. 
            plot(ax,[min(delay(i,:)):1:-100+zeroDelay 100+zeroDelay:1:max(delay(i,:))],...
                result(abs([min(delay(i,:)):1:-100+zeroDelay 100+zeroDelay:1:max(delay(i,:))])),'r','DisplayName','');   
        case 'stretched-exp'
            g = fittype(@(A,B,beta,x) A*exp(-B*abs(x).^beta) );                   
            B0 = 1/max(x);
            beta0 = 1;
            A0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay)));
            [result,gof,~] = fit(x',ys',g,'StartPoint',[A0 B0 beta0]); 
            fitTau(i) = result.beta;
            [se]= getStandardErrorsFromFit(result,gof,'method1'); 
            B = result.B;
            beta = result.beta;
            pa1(i) = result.A;
            pa2(i) = result.B;
            pa3(i) = result.beta;
            pa1Err(i) = se(1);
            pa2Err(i) = se(2);
            pa3Err(i) = se(3);
            fitPeak(i) = result.A; 
            fitPeakErr(i) = se(1); 
            
            % from the paper and supplement Caputo, D., et al. (2018).
            %Topological order and thermal equilibrium in polariton condensates. 
            %Nature Materials, 17(2), 145?151. https://doi.org/10.1038/NMAT5039
            % set 2x = t in Eq. 1 and 2
%             meanTime = B^(-1/beta);
%             fitTau(i) = meanTime/beta * gamma(1/beta);
            
            % from the definition of coherence time in https://www.rp-photonics.com/coherence_time.html
            % coherence time is given by integrating over the square of
            % g^(1)(t) from -inf to +inf. This can be done by wolfram alpha
            % and gives the result:
            fitTau(i) = 2^((-1 + beta)/beta)*B^(-1/beta)*gamma(1+ 1/beta);
            % get error from monte carlo error propagation
            [fitTauRand] = deal(zeros(1000,1));           
            for u = 1:1000
                Brand = normrnd(B,se(2));
                betarand = normrnd(beta,se(3));
                fitTauRand(u) = 2^((-1 + betarand)/betarand)*Brand^(-1/betarand)*gamma(1+ 1/betarand);
            end
            tauErr(i) = std(fitTauRand);            
            plot(ax,min(delay(i,:)):1:max(delay(i,:)),result(abs(min(delay(i,:)):1:max(delay(i,:)))),'r','DisplayName','');           
         case 'lorentz'
             % starting values 
            p02 = zeroDelay;
            p03 = 0.5*max(x);
            C = mean(ys(end-5:end)); % saturation value 
            peakHeight = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay))); 
            p01 = (peakHeight- C)*p03;
            p0 = [p01 p02 p03 C];
             [result,params,RESNORM, RESIDUAL, JACOBIAN] = lorentzfit(x',ys',p0);
            p1 = params(1);  p2 = params(2);  p3 = params(3); C = params(4);
            pa1(i) = p1;
            pa2(i) = p2;
            pa3(i) = p3;
            pa4(i) = p4;
            fitTau(i) = 2*sqrt(p3);  %FWHM
            fitPeak(i) = p1/p3 + C;  
            plot(ax,delay(i,:),result,'r','DisplayName','');                   
        case 'voigt' 
%                     fo = fitoptions('Method','NonlinearLeastSquares', ...
%                             'StartPoint',[50, 5]); %gamma, sigma
%                     zeroDelay = 200;  %set this...  
%                     ysFit = ys/max(ys);
%                     ft = fittype( @(gam,sig,peakVal,x) myvoigt(x, peakVal, gam, sig) ,'problem','peakVal','options',fo);    
%                     [res,gof,~] = fit(x',ysFit,ft,'problem',zeroDelay);
%                     resPlot= res(min(delay(i,:)):1:max(delay(i,:)));
%                     resPlot = resPlot * max(ys);
%                     plot(ax,min(delay(i,:)):1:max(delay(i,:)),resPlot,'r','DisplayName','');              
            initGuess1 = [0,30, 1]; %peak x, gamma, sigma
            [estimates1, model1] = voigtfit(x, ys, initGuess1, [min(x), max(x)]);
            gam = estimates1(2);
            sigma = estimates1(3);
            result = myvoigt(x, 0, gam, sigma );%estimates1(1)
            result = result';
            c = abs(result\ys');
            result = result*c;             
            plot(ax,x,result,'r','DisplayName','');
            FWHM_Gauss = 2*sigma*sqrt(2*log(2));
            FWHM_Lorentz = 2*gam;         
            pa1(i) = estimates1(i);
            pa2(i) = gam;
            pa3(i) = sigma; 
            fitTau(i) = 0.5346*FWHM_Lorentz + sqrt(0.2166*FWHM_Lorentz^2 + FWHM_Gauss^2);
            fitPeak(i) = max(result);
            % https://en.wikipedia.org/wiki/Voigt_profile
            % Olivero, J. J.; R. L. Longbothum (February 1977). 
%                     "Empirical fits to the Voigt line width: A brief review". 
%                     Journal of Quantitative Spectroscopy and Radiative Transfer. 
%                     17 (2): 233???236. Bibcode:1977JQSRT..17..233O. doi:10.1016/0022-4073(77)90161-3. ISSN 0022-4073.
        case 'noFit'
            fitPeak(i) = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay)));
    end %fitType
    
    if ~(strcmp(fitType,'noFit') || strcmp(fitType,'lorentz') || strcmp(fitType,'voigt'))
        sse(i) = gof.sse;
        rsquare(i) = gof.rsquare;
        adjrsquare(i) = gof.adjrsquare;
        rmse(i) = gof.rmse;
    end
    
end %iloop
hold off;
if strcmp(fitType,'noFit')
    legend('location','southeast','Fontsize',20);
else
    f=get(ax,'Children');
    index = length(f)-((1:length(I))-1).*2;
    legend(f(index),'location','southeast','Fontsize',20);
end
xlabel(ax,['Delay (' xUnit ')']);
if logplot
    ax.YScale = 'log';
end
fig = figure(1);              
set(fig,'Color','w','Units','centimeters','Position',[1,1,45,30],'PaperPositionMode','auto');
graphicsSettings;
set(ax,'FontSize',38);
%% Write figure to file and close it
if ~isempty(filename)
    savefig(fig,filename);
    print([filename '.png'],'-dpng');
    close all;
end

%% plot fit results vs ring radius and thickness
if any(fitTau)
    YrPlot = Yr(:,1);
    [YrPlot,Ir]= sort(YrPlot);
    fitTau = real(fitTau(Ir));
    tauErr = tauErr(Ir);
    errorbar(YrPlot,fitTau,tauErr,'o-','Linewidth',2);
    xlabel('r_{ps} set for postselection');
    switch fitType
        case 'gauss'
            ylabel(['\tau_c of Gaussian (' xUnit ')']);
        case 'voigt'
            ylabel(['FWHM of Voigt Profile (' xUnit ')']);
        case 'exponential'
            ylabel(['\tau_c of exp. function (' xUnit ')']);
        case 'exp2'
            ylabel(['\tau_c of exp. function (' xUnit ')']);
        case 'exp3'
            ylabel(['\tau_c of exp. function (' xUnit ')']);
        case 'sech'
            ylabel(['\tau_c of sech function (' xUnit ')']);
         case 'lorentz'
            ylabel(['FWHM of Lorentz Profile (' xUnit ')']);
         case 'power-law'
            ylabel('\alpha of power-law Profile');
        case 'stretched-exp'
            ylabel('\tau_c of stretched-exp Profile');
    end
    graphicsSettings;
    ax = gca;
    set(ax,'FontSize',30);
    title(typestr);
    savefig([filename '-fitWidthsVsRadius.fig']);
    print([filename '-fitWidthsVsRadius.png'],'-dpng');
    close all;
end

if any(fitPeak)
    YrPlot = Yr(:,1);
    [YrPlot,Ir]= sort(YrPlot);
    fitPeak = real(fitPeak(Ir));
    plot(YrPlot,fitPeak,'o-');
    xlabel('r_{ps} set for postselection');
    switch fitType
        case 'power-law'
            ylabel('A of power-law Profile');
        otherwise
            ylabel('Peakheight');
    end
    graphicsSettings;
    ax = gca;
    set(ax,'FontSize',30);
    title(typestr);
    savefig([filename '-PeaksVsRadius.fig']);
    print([filename '-PeaksVsRadius.png'],'-dpng');
    close all;
end

%% save fitresults
save([filename '-fitresults.mat'],'YrPlot','fitTau','tauErr','fitPeak','fitPeakErr',...
    'sse','rsquare','adjrsquare','rmse','pa1','pa2','pa3','pa4','pa1Err','pa2Err','pa3Err','pa4Err');
T = table(YrPlot,fitTau,tauErr,fitPeak,fitPeakErr,sse,rsquare,adjrsquare,rmse,pa1,pa2,pa3,pa4,pa1Err,pa2Err,pa3Err,pa4Err);
writetable(T,[filename '-fitresults.txt'],'WriteVariableNames',true);

% 
% if any(fitTau)
%     YtPlot = Yt(:,1);
%     [YtPlot,It]= sort(YrPlot);
%     fitTau = fitTau(It);
%     tauError = tauError(It);
%     errorbar(YtPlot,fitTau,tauError,'o-','Linewidth',2);
%     xlabel('Ring Thickness set for postselection');
%     switch fitType
%         case 'gauss'
%             ylabel(['\tau_c of Gaussian (' xUnit ')']);
%         case 'voigt'
%             ylabel(['FWHM of Voigt Profile (' xUnit ')']);
%          case 'exponential'
%             ylabel(['\tau of exp. function (' xUnit ')']);
%     end
%     graphicsSettings;
%     title([filename '-' typestr '-' fitType '-fitWidthsVsThickness']);
%     savefig([filename '-' typestr '-' fitType '-fitWidthsVsThickness.fig']);
%     print([filename '-' typestr '-' fitType '-fitWidthsVsThickness.png'],'-dpng');
%     close all;
% end
% 
% if any(fitPeak)
%     fitPeak = fitPeak(It);
%     plot(YtPlot,fitPeak,'o-');
%     xlabel('Ring Thickness set for postselection');
%     ylabel('Peakheight');
%     graphicsSettings;
%     title([filename '-' typestr '-' fitType '-PeaksVsThickness']);
%     savefig([filename '-' typestr '-' fitType '-PeaksVsThickness.fig']);
%     print([filename '-' typestr '-' fitType '-PeaksVsThickness.png'],'-dpng');
%     close all;
% end
% plot(Yt(:,1),fitWidth,'o-');
% plot(Yt(:,1),fitPeak,'o-');
    end %plotStuff

end

