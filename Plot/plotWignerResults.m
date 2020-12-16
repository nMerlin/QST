function plotWignerResults(listOfParams,varargin)
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
defaultSmooth = false;
addParameter(p,'Smooth',defaultSmooth,@islogical);
defaultPlotrelative = false; %plots the data relative to the target photon number
addParameter(p,'Plotrelative',defaultPlotrelative,@islogical);
parse(p,varargin{:});
c = struct2cell(p.Results);
[fitType,plotrelative,range,remMod,smooth,varyAPS,xUnit,zeroDelay] = c{:};



%% Create folder 'figures-fig'
if ~exist([pwd 'Wigner-figures-fig'],'dir')
    mkdir('Wigner-figures-fig')
end

% Constants
figurepath = 'Wigner-figures-fig/';

load('Photonnumbers.mat','nPsFast','nPsSlow','nTg'); %loads the photon numbers of each channel without postselection,...
%which was created with plotSeriesPostselections
%% Gather data
[delay,Yr,Yt,Qs,Ps,varQs,varPs,discN,meanPh,meanR,varPh,varR,meanAbsPh] = deal([]);
for iParams = 1:length(listOfParams)
    selParams = listOfParams(iParams);
    selStr = selParamsToStr(selParams);
    foldername = ['Wignerplots-',selStr,'-remMod-',...
            num2str(remMod),'-range-',num2str(range),'-varyAPS-',num2str(varyAPS)];
    load([foldername '\Wignerresults-smooth-' num2str(smooth) '.mat'],'Delay','Q','P','varQ','varP','n',...
        'meanPhases','meanAmps','varPhases','varAmps','meanAbsPhases');
    delay(iParams,:) = Delay; 
    H = length(delay);
    Radii = ones(H,1) * selParams.Position(1);
    Thicknesses = ones(H,1) * selParams.Position(2);
    Yr(iParams,:) = Radii;
    Yt(iParams,:) = Thicknesses;
    Qs(iParams,:) = Q;
    Ps(iParams,:) = P;
    varQs(iParams,:) = varQ;
    varPs(iParams,:) = varP;
    discN(iParams,:) = n;
    meanPh(iParams,:) = meanPhases;
    meanR(iParams,:) = meanAmps; 
    varPh(iParams,:) = varPhases;
    varR(iParams,:) = varAmps;
    meanAbsPh(iParams,:) = meanAbsPhases;
end
[~,I] = sort(delay(:,1)); % Sort for Radii
delay = delay(I,:); %delays
Yr = Yr(I,:); %Radii
Yt = Yt(I,:); %thicknesses
Qs= Qs(I,:);
Ps= Ps(I,:);
varQs = varQs(I,:);
varPs = varPs(I,:);
discN = discN(I,:);
meanPh = meanPh(I,:);
meanR = meanR(I,:);varPh = varPh(I,:);varR = varR(I,:);meanAbsPh = meanAbsPh(I,:);
[fitTau,fitPeak,tauError] = deal(zeros(length(I),1));

%% Plot
typestrVector = {'R','Q','P','RLog','meanPh','meanAbsPh','discN','varR','varPh','varQ','varP'};
%typestrVector = {'R','discN','varR'};
%typestrVector = {'R'};
for typeI = 1:length(typestrVector)
    plotStuff(cell2mat(typestrVector(typeI)))
end

    function [] = plotStuff(typestr)
%% Create figure
fig = figure;
filename = [figurepath,typestr,'-',fitType,'-remMod-',...
            num2str(remMod),'-range-',num2str(range),'-varyAPS-',num2str(varyAPS),'-smooth-',num2str(smooth) '.fig'];
%formatFigA5(fig);
switch typestr
    case 'R'
        figure(1);
        ax = gca;
        for i = 1:length(I)
            plot(ax,delay(i,:),meanR(i,:),'o-','DisplayName',['r = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
            hold on;
            x = delay(i,:);
            if ~any(~isnan(meanR(i,:))) %if there are only nans
                meanR(i,:) = zeros(1,H);
            end
%             ys = transpose(csaps(x,meanR(i,:),0.00001,x));  
%             ys = ys';
            ys = meanR(i,:);
            switch fitType
                case 'gauss'
                    g = fittype(@(a,b,c,d,x) a*exp(-pi/2*((x-b)/c).^2)+d);                   
                    b0 = zeroDelay;
                    c0 = 0.5*max(x);
                    d0 = mean(ys(end-5:end)); % saturation value
                    a0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay))) - d0;
                    [res,gof,~] = fit(x',ys',g,'StartPoint',[a0 b0 c0 d0]); 
                    fitTau(i) = res.c;
                    level = 2*tcdf(-1,gof.dfe);
                    m = confint(res,1-level); 
                    tauError(i) = m(end,3) - res.c;
                    fitPeak(i) = res(res.b);  
                    plot(ax,min(delay(i,:)):1:max(delay(i,:)),res(min(delay(i,:)):1:max(delay(i,:))),'r','DisplayName','');
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
                    [res,gof,~] = fit(abs(x)',log(ys)',ex,'StartPoint',[m0 c0]); 
                    fitTau(i) = res.c;
                    level = 2*tcdf(-1,gof.dfe);
                    m = confint(res,1-level); 
                    tauError(i) = m(end,2) - res.c;         
                    fitPeak(i) = sig*exp(res(b0))+d;  
                    plot(ax,min(delay(i,:)):1:max(delay(i,:)),sig*exp(res(abs(min(delay(i,:)):1:max(delay(i,:)))))+d,'r','DisplayName',''); 
                case 'sech'
                    se = fittype(@(a,b,c,d,x) a./(cosh(-pi/2*(x-b)/c)).^2+d);                   
                    b0 = zeroDelay;
                    c0 = 0.5*max(x);
                    d0 = mean(ys(end-5:end)); % saturation value
                    a0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay))) - d0;
                    [res,gof,~] = fit(x',ys',se,'StartPoint',[a0 b0 c0 d0]); 
                    fitTau(i) = res.c;
                    level = 2*tcdf(-1,gof.dfe);
                    m = confint(res,1-level); 
                    tauError(i) = m(end,3) - res.c;
                    fitPeak(i) = res(res.b);  
                    plot(ax,min(delay(i,:)):1:max(delay(i,:)),res(min(delay(i,:)):1:max(delay(i,:))),'r','DisplayName','');
                 case 'lorentz'
                     % starting values 
                    p02 = zeroDelay;
                    p03 = 0.5*max(x);
                    C = mean(ys(end-5:end)); % saturation value 
                    peakHeight = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay))); 
                    p01 = (peakHeight- C)*p03;
                    p0 = [p01 p02 p03 C];
                     [res,params,RESNORM, RESIDUAL, JACOBIAN] = lorentzfit(x',ys',p0);
                    p1 = params(1);  p2 = params(2);  p3 = params(3); C = params(4);
                    fitTau(i) = 2*sqrt(p3);  %FWHM
                    fitPeak(i) = p1/p3 + C;  
                    plot(ax,delay(i,:),res,'r','DisplayName','');                   
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
                    gamma = estimates1(2);
                    sigma = estimates1(3);
                    res = myvoigt(x, 0, gamma, sigma );%estimates1(1)
                    res = res';
                    c = abs(res\ys');
                    res = res*c;             
                    plot(ax,x,res,'r','DisplayName','');
                    FWHM_Gauss = 2*sigma*sqrt(2*log(2));
                    FWHM_Lorentz = 2*gamma;
                    fitTau(i) = 0.5346*FWHM_Lorentz + sqrt(0.2166*FWHM_Lorentz^2 + FWHM_Gauss^2);
                    fitPeak(i) = max(res);
                    % https://en.wikipedia.org/wiki/Voigt_profile
                    % Olivero, J. J.; R. L. Longbothum (February 1977). 
%                     "Empirical fits to the Voigt line width: A brief review". 
%                     Journal of Quantitative Spectroscopy and Radiative Transfer. 
%                     17 (2): 233�236. Bibcode:1977JQSRT..17..233O. doi:10.1016/0022-4073(77)90161-3. ISSN 0022-4073.
                case 'noFit'
                    Index = find(delay(i,:)>= -30 & delay(i,:)<= 30);
                    fitPeak(i) = mean(meanR(i,Index-2:Index+2)); 
            end              
        end
        hold off;
        if strcmp(fitType,'noFit')
            legend('location','best');
        else
            f=get(ax,'Children');
            index = length(f)-((1:length(I))-1).*2;
            legend(f(index),'location','best');
        end
        xlabel(ax,['Delay (' xUnit ')']); 
        ylabel(ax,'mean amplitude <r>'); 
        fig = figure(1);       
    case 'RLog'
        figure(1);
        ax = gca;
        for i = 1:length(I)
            semilogy(ax,delay(i,:),meanR(i,:),'o-','DisplayName',['r = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
            hold on;   
        end
        hold off;
        legend('location','best');
        xlabel(ax,['Delay (' xUnit ')']); 
        ylabel(ax,'mean amplitude <r>'); 
        fig = figure(1);   
    case 'Q'
        figure(1);
        ax = gca;
        for i = 1:length(I)
            plot(ax,delay(i,:),Qs(i,:),'o-','DisplayName',['r = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
            hold on;   
        end
        hold off;
        legend('location','best');
        xlabel(ax,['Delay (' xUnit ')']); 
        ylabel(ax,'<Q>'); 
        fig = figure(1);              
    case 'P'
        figure(1);
        ax = gca;
        for i = 1:length(I)
            plot(ax,delay(i,:),Ps(i,:),'o-','DisplayName',['r = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
            hold on;   
        end
        hold off;
        legend('location','best');
        xlabel(ax,['Delay (' xUnit ')']); 
        ylabel(ax,'<P>'); 
        fig = figure(1);
    case 'meanPh'
        figure(1);
        ax = gca;
        for i = 1:length(I)
            plot(ax,delay(i,:),meanPh(i,:),'o-','DisplayName',['r = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
            hold on;   
        end
        hold off;
        legend('location','best');
        xlabel(ax,['Delay (' xUnit ')']); 
        ylabel(ax,'<\phi> (rad)'); 
        fig = figure(1);
    case 'meanAbsPh'
        figure(1);
        ax = gca;
        for i = 1:length(I)
            plot(ax,delay(i,:),meanAbsPh(i,:),'o-','DisplayName',['r = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
            hold on;   
        end
        hold off;
        legend('location','best');
        xlabel(ax,['Delay (' xUnit ')']); 
        ylabel(ax,'<|\phi|> (rad)'); 
        fig = figure(1);
     case 'varPh'
         figure(1);
         ax = gca;
        for i = 1:length(I)
            plot(ax,delay(i,:),varPh(i,:),'o-','DisplayName',['r = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
            hold on;
            Index = find(delay(i,:)>= -30 & delay(i,:)<= 30);
            fitPeak(i) = mean(varPh(i,Index)); 
        end
        hold off;
        legend('location','southeast');
        ylabel(ax,'Variance(\phi) (rad)');
        xlabel(ax,['Delay (' xUnit ')']);
        fig = figure(1);
    case 'varR'
        figure(1);
        ax= gca;
        for i = 1:length(I)
            plot(delay(i,:),varR(i,:),'o-','DisplayName',['r = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
            hold on;
            x = delay(i,:);
            if ~any(~isnan(varR(i,:))) %if there are only nans
                varR(i,:) = zeros(1,H);
            end
%             ys = transpose(csaps(x,varR(i,:),0.00001,x));  
%             ys = ys';
            ys = varR(i,:);
            switch fitType
                case 'gauss'
                   g = fittype(@(a,b,c,d,x) a*exp(-pi/2*((x-b)/c).^2)+d);                   
                    b0 = zeroDelay;
                    c0 = 0.5*max(x);
                    d0 = mean(ys(end-5:end)); % saturation value
                    a0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay))) - d0;
                    [res,gof,~] = fit(x',ys',g,'StartPoint',[a0 b0 c0 d0]); 
                    fitTau(i) = res.c;
                    level = 2*tcdf(-1,gof.dfe);
                    m = confint(res,1-level); 
                    tauError(i) = m(end,3) - res.c;
                    fitPeak(i) = res(res.b);  
                    plot(ax,min(delay(i,:)):1:max(delay(i,:)),res(min(delay(i,:)):1:max(delay(i,:))),'r','DisplayName','');
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
                    [res,gof,~] = fit(abs(x)',log(ys)',ex,'StartPoint',[m0 c0]); 
                    fitTau(i) = res.c;
                    level = 2*tcdf(-1,gof.dfe);
                    m = confint(res,1-level); 
                    tauError(i) = m(end,2) - res.c;         
                    fitPeak(i) = sig*exp(res(b0))+d;  
                    plot(ax,min(delay(i,:)):1:max(delay(i,:)),sig*exp(res(abs(min(delay(i,:)):1:max(delay(i,:)))))+d,'r','DisplayName','');  
                case 'sech'
                    se = fittype(@(a,b,c,d,x) a./(cosh(-pi/2*(x-b)/c)).^2+d);                   
                    b0 = zeroDelay;
                    c0 = 0.5*max(x);
                    d0 = mean(ys(end-5:end)); % saturation value
                    a0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay))) - d0;
                    [res,gof,~] = fit(x',ys',se,'StartPoint',[a0 b0 c0 d0]); 
                    fitTau(i) = res.c;
                    level = 2*tcdf(-1,gof.dfe);
                    m = confint(res,1-level); 
                    tauError(i) = m(end,3) - res.c;
                    fitPeak(i) = res(res.b);  
                    plot(ax,min(delay(i,:)):1:max(delay(i,:)),res(min(delay(i,:)):1:max(delay(i,:))),'r','DisplayName','');
                 case 'lorentz'
                      % starting values 
                    p02 = zeroDelay;
                    p03 = 0.5*max(x);
                    C = mean(ys(end-5:end)); % saturation value 
                    peakHeight = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay))); 
                    p01 = (peakHeight- C)*p03;
                    p0 = [p01 p02 p03 C];
                     [res,params,RESNORM, RESIDUAL, JACOBIAN] = lorentzfit(x',ys',p0);
                    p1 = params(1);  p2 = params(2);  p3 = params(3); C = params(4);
                    fitTau(i) = 2*sqrt(p3);  %FWHM
                    fitPeak(i) = p1/p3 + C;  
                    plot(ax,delay(i,:),res,'r','DisplayName','');                   
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
                    gamma = estimates1(2);
                    sigma = estimates1(3);
                    res = myvoigt(x, 0, gamma, sigma );%estimates1(1)
                    res = res';
                    c = abs(res\ys');
                    res = res*c;             
                    plot(ax,x,res,'r','DisplayName','');
                    FWHM_Gauss = 2*sigma*sqrt(2*log(2));
                    FWHM_Lorentz = 2*gamma;
                    fitTau(i) = 0.5346*FWHM_Lorentz + sqrt(0.2166*FWHM_Lorentz^2 + FWHM_Gauss^2);
                    fitPeak(i) = max(res);
                    % https://en.wikipedia.org/wiki/Voigt_profile
                    % Olivero, J. J.; R. L. Longbothum (February 1977). 
%                     "Empirical fits to the Voigt line width: A brief review". 
%                     Journal of Quantitative Spectroscopy and Radiative Transfer. 
%                     17 (2): 233�236. Bibcode:1977JQSRT..17..233O. doi:10.1016/0022-4073(77)90161-3. ISSN 0022-4073.
                case 'noFit'
                    Index = find(delay(i,:)>= -30 & delay(i,:)<= 30);
                    fitPeak(i) = mean(meanR(i,Index-2:Index+2)); 
            end              
        end
        hold off;
        if strcmp(fitType,'noFit')
            legend('location','best');
        else
            f=get(ax,'Children');
            index = length(f)-((1:length(I))-1).*2;
            legend(f(index),'location','best');
        end
%             Index = find(delay(i,:)>= -30 & delay(i,:)<= 30);
%             fitPeak(i) = mean(varR(i,Index)); 
        ylabel('Variance(r)');
        xlabel(['Delay (' xUnit ')']);
        fig = figure(1);
     case 'varQ'
         figure(1);
         ax = gca;
        for i = 1:length(I)
            plot(delay(i,:),varQs(i,:),'o-','DisplayName',['r = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
            hold on;
            Index = find(delay(i,:)>= -30 & delay(i,:)<= 30);
            fitPeak(i) = mean(varQs(i,Index)); 
        end;
        hold off;
        legend('location','southeast');
        ylabel('Variance in Q');
        xlabel(['Delay (' xUnit ')']);
        fig = figure(1);
     case 'varP'
         figure(1);
         ax = gca;
        for i = 1:length(I)
            plot(delay(i,:),varPs(i,:),'o-','DisplayName',['r = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
            hold on;
            Index = find(delay(i,:)>= -30 & delay(i,:)<= 30);
            fitPeak(i) = mean(varPs(i,Index)); 
        end;
        hold off;
        legend('location','southeast');
        ylabel('Variance in P');
        xlabel(['Delay (' xUnit ')']);
        fig = figure(1);
    case 'discN'
        figure(1);
        ax = gca;
        for i = 1:length(I)
            plot(ax,delay(i,:),discN(i,:),'o-','DisplayName',['r = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
            hold on;
             x = delay(i,:);
            if ~any(~isnan(discN(i,:))) %if there are only nans
                discN(i,:) = zeros(1,H);
            end
%             ys = transpose(csaps(x,meanR(i,:),0.00001,x));  
%             ys = ys';
            ys = discN(i,:);
            switch fitType
                case 'gauss'
                    g = fittype(@(a,b,c,d,x) a*exp(-pi/2*((x-b)/c).^2)+d);                   
                    b0 = zeroDelay;
                    c0 = 0.5*max(x);
                    d0 = mean(ys(end-5:end)); % saturation value
                    a0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay))) - d0;
                    [res,gof,~] = fit(x',ys',g,'StartPoint',[a0 b0 c0 d0]); 
                    fitTau(i) = res.c;
                    level = 2*tcdf(-1,gof.dfe);
                    m = confint(res,1-level); 
                    tauError(i) = m(end,3) - res.c;
                    fitPeak(i) = res(res.b);  
                    plot(ax,min(delay(i,:)):1:max(delay(i,:)),res(min(delay(i,:)):1:max(delay(i,:))),'r','DisplayName','');
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
                    [res,gof,~] = fit(abs(x)',log(ys)',ex,'StartPoint',[m0 c0]); 
                    fitTau(i) = res.c;
                    level = 2*tcdf(-1,gof.dfe);
                    m = confint(res,1-level); 
                    tauError(i) = m(end,2) - res.c;         
                    fitPeak(i) = sig*exp(res(b0))+d;  
                    plot(ax,min(delay(i,:)):1:max(delay(i,:)),sig*exp(res(abs(min(delay(i,:)):1:max(delay(i,:)))))+d,'r','DisplayName','');
                case 'sech'
                    se = fittype(@(a,b,c,d,x) a./(cosh(-pi/2*(x-b)/c)).^2+d);                   
                    b0 = zeroDelay;
                    c0 = 0.5*max(x);
                    d0 = mean(ys(end-5:end)); % saturation value
                    a0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay))) - d0;
                    [res,gof,~] = fit(x',ys',se,'StartPoint',[a0 b0 c0 d0]); 
                    fitTau(i) = res.c;
                    level = 2*tcdf(-1,gof.dfe);
                    m = confint(res,1-level); 
                    tauError(i) = m(end,3) - res.c;
                    fitPeak(i) = res(res.b);  
                    plot(ax,min(delay(i,:)):1:max(delay(i,:)),res(min(delay(i,:)):1:max(delay(i,:))),'r','DisplayName','');
                 case 'lorentz' 
                     % starting values 
                    p02 = zeroDelay;
                    p03 = 0.5*max(x);
                    C = mean(ys(end-5:end)); % saturation value 
                    peakHeight = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay))); 
                    p01 = (peakHeight- C)*p03;
                    p0 = [p01 p02 p03 C];
                     [res,params,RESNORM, RESIDUAL, JACOBIAN] = lorentzfit(x',ys',p0);
                    p1 = params(1);  p2 = params(2);  p3 = params(3); C = params(4);
                    fitTau(i) = 2*sqrt(p3);  %FWHM
                    fitPeak(i) = p1/p3 + C; 
                    plot(ax,delay(i,:),res,'r','DisplayName','');                   
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
                    gamma = estimates1(2);
                    sigma = estimates1(3);
                    res = myvoigt(x, 0, gamma, sigma );%estimates1(1)
                    res = res';
                    c = abs(res\ys');
                    res = res*c;             
                    plot(ax,x,res,'r','DisplayName','');
                    FWHM_Gauss = 2*sigma*sqrt(2*log(2));
                    FWHM_Lorentz = 2*gamma;
                    fitTau(i) = 0.5346*FWHM_Lorentz + sqrt(0.2166*FWHM_Lorentz^2 + FWHM_Gauss^2);
                    fitPeak(i) = max(res);
                    % https://en.wikipedia.org/wiki/Voigt_profile
                    % Olivero, J. J.; R. L. Longbothum (February 1977). 
%                     "Empirical fits to the Voigt line width: A brief review". 
%                     Journal of Quantitative Spectroscopy and Radiative Transfer. 
%                     17 (2): 233�236. Bibcode:1977JQSRT..17..233O. doi:10.1016/0022-4073(77)90161-3. ISSN 0022-4073.
                case 'noFit'
                    Index = find(delay(i,:)>= -30 & delay(i,:)<= 30);
                    fitPeak(i) = mean(meanR(i,Index-2:Index+2)); 
            end        
        end
        hold off;
        if strcmp(fitType,'noFit')
            legend('location','best');
        else
            f=get(ax,'Children');
            index = length(f)-((1:length(I))-1).*2;
            legend(f(index),'location','best');
        end
        xlabel(['Delay (' xUnit ')']);
        ylabel('Photon number');
        fig = figure(1);
end
set(fig,'Color','w','Units','centimeters','Position',[1,1,45,30],'PaperPositionMode','auto');
graphicsSettings;
set(ax,'FontSize',38);
%% Write figure to file and close it
if ~isempty(filename)
    savefig(fig,filename);
    print([filename '.png'],'-dpng');
    ax.YScale = 'log';
    savefig(fig,[filename '-log.fig']);
    print([filename 'log.png'],'-dpng');
    close all;
end

% %% plot fit results vs ring radius and thickness
if any(fitTau)
    YrPlot = Yr(:,1);
    [YrPlot,Ir]= sort(YrPlot);
    fitTau = real(fitTau(Ir));
    tauError = tauError(Ir);
    errorbar(YrPlot,fitTau,tauError,'o-','Linewidth',2);
    xlabel('A_{ps} set for postselection');
    switch fitType
        case 'gauss'
            ylabel(['\tau_c of Gaussian (' xUnit ')']);
        case 'voigt'
            ylabel(['FWHM of Voigt Profile (' xUnit ')']);
        case 'exponential'
            ylabel(['\tau_c of exp. function (' xUnit ')']);
        case 'sech'
            ylabel(['\tau_c of sech function (' xUnit ')']);
         case 'lorentz'
            ylabel(['FWHM of Lorentz Profile (' xUnit ')']);
    end
    graphicsSettings;
    title(typestr);
    savefig([filename '-fitWidthsVsRadius.fig']);
    print([filename '-fitWidthsVsRadius.png'],'-dpng');
    close all;
end

if any(fitPeak)
    fitPeak = real(fitPeak(Ir));
    plot(YrPlot,fitPeak,'o-');
    xlabel('A_{ps} set for postselection');
    ylabel('Peakheight');
    graphicsSettings;
    title(typestr);
    savefig([filename '-PeaksVsRadius.fig']);
    print([filename '-PeaksVsRadius.png'],'-dpng');
    close all;
end
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
    end

end

