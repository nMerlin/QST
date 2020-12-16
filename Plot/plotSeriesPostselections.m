function plotSeriesPostselections(listOfParams,varargin)
%PLOTSERIESPOSTSELECTIONS

%% Validate and parse input arguments
p = inputParser;
defaultFilename = '';
addParameter(p,'Filename',defaultFilename,@isstr);
defaultType = 'Amplitude';
addParameter(p,'Type',defaultType,@isstr);
defaultXUnit = 'fs';
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
parse(p,varargin{:});
c = struct2cell(p.Results);
[filename,fitType,plotrelative,range,remMod,typestr,varyAPS,xUnit,zeroDelay] = c{:};
filename = [filename '-varyAPS-' num2str(varyAPS) '-remMod-' num2str(remMod) '-fitType-' fitType '-PlotRelative-' num2str(plotrelative) '.fig'];

%% Create folder 'figures-fig'
if ~exist([pwd 'figures-fig'],'dir')
    mkdir('figures-fig')
end

% Constants
figurepath = 'figures-fig/';

%% Gather data
[delay,Yr,Yt,discAmpl,discMeanVar,varQ,varP,discN,g2vals,g2std,nTg,nPsFast,nPsSlow] = deal([]);
sigmas = zeros(length(listOfParams),1);
sigmaConf = zeros(length(listOfParams),2);
for iParams = 1:length(listOfParams)
    selParams = listOfParams(iParams);
    
    % From tables
    A = seriesRead3ChTable(selParams,'VaryAPS',varyAPS,'RemoveModulation',remMod,'Range',range);
    H = height(A);
    Radii = ones(H,1) * selParams.Position(1);
    Thicknesses = ones(H,1) * selParams.Position(2);
    delay(iParams,:) = A.Delay;
    [delay(iParams,:),I] = sort(delay(iParams,:)); % Sort for Delays
    Yr(iParams,:) = Radii;
    Yt(iParams,:) = Thicknesses;
    discAmpl(iParams,:) = A.discAmpl;
    discAmpl(iParams,:) = discAmpl(iParams,I);
    discMeanVar(iParams,:) = A.discMeanVar;
    discMeanVar(iParams,:) = discMeanVar(iParams,I);
    varQ(iParams,:) = A.varQ;
    varQ(iParams,:) = varQ(iParams,I);
    varP(iParams,:) = A.varP;
    varP(iParams,:) = varP(iParams,I);
    discN(iParams,:) = A.discN;
    discN(iParams,:) = discN(iParams,I);
    g2vals(iParams,:) = A.g2;
    g2vals(iParams,:) = g2vals(iParams,I);
    g2std(iParams,:) = A.g2std;
    g2std(iParams,:) = g2std(iParams,I);
    g2std(iParams,:) = A.g2std;
    g2std(iParams,:) = g2std(iParams,I);
    nTg(iParams,:) = A.nTg;
    nTg(iParams,:) = nTg(iParams,I);
    nPsFast(iParams,:) = A.nPsFast;
    nPsFast(iParams,:) = nPsFast(iParams,I);
    nPsSlow(iParams,:) = A.nPsSlow;
    nPsSlow(iParams,:) = nPsSlow(iParams,I);
    
%     nX1(iParams,:) = A.nX1;
%     nX1(iParams,:) = nX1(iParams,I);
%     nX2(iParams,:) = A.nX2;
%     nX2(iParams,:) = nX2(iParams,I);
%     nX3(iParams,:) = A.nX3;
%     nX3(iParams,:) = nX3(iParams,I);
    
    % From DelayMeanVarX Plots
    if strcmp(typestr,'MeanVarSigma') || ...
            strcmp(typestr,'ThicknessMeanVarSigma')
        selstr = selParamsToStr(selParams);
        filelist = dir([figurepath,'*-DelayMeanVarX-',selstr,'.fig']);
        filelist = {filelist.name};
        fig = openfig([figurepath,filelist{1}]);
        figData = get(gca,'Children');
        fitStr = strjoin(figData(1).String);
        toks = regexpi(fitStr,'s =\s*([\d.-]*)','tokens');
        sigmas(iParams) = str2double(cell2mat(toks{1}));
        toks = regexpi(fitStr,'s =[\d.\s]*\(([\d.-]*)','tokens');
        sigmaConf(iParams,1) = str2double(cell2mat(toks{1}));
        toks = regexpi(fitStr,'s =[\d.\s]*\([\d.]*,\s([\d.-]*)','tokens');
        sigmaConf(iParams,2) = str2double(cell2mat(toks{1}));
        close(fig);
    end
    
    % From DelayDiscAmpl Plots
end
[~,I] = sort(delay(:,1)); % Sort for Radii
delay = delay(I,:); %delays
c = 299792458; % in m/s
delay = delay - 2*zeroDelay/1000/c*10^12;
Yr = Yr(I,:); %Radii
Yt = Yt(I,:); %thicknesses
discAmpl = discAmpl(I,:);
discMeanVar = discMeanVar(I,:);
discN = discN(I,:);
g2vals = g2vals(I,:);
g2std = g2std(I,:);
nTg = nTg(1,:);
nPsFast = nPsFast(1,:);
nPsSlow = nPsSlow(1,:);
meanNPsFast = mean(nPsFast);
meanNPsSlow = mean(nPsSlow);
meanNTg = mean(nTg);
save('Photonnumbers.mat','nPsFast','nPsSlow','nTg','meanNPsFast','meanNPsSlow','meanNTg');
[fitTau,fitPeak,tauError] = deal(zeros(length(I),1));

if plotrelative
    discAmpl = discAmpl./sqrt(nTg);
    discN = discN./nTg;
    discMeanVar = discMeanVar./nTg;
    varQ = varQ./nTg;
    varP = varP./nTg;
    sel = 'k'; %the selection radius or variable
else
    sel = 'r';
end

%% Create figure
fig = figure;
%formatFigA5(fig);
switch typestr
    case 'Amplitude'
        figure(1);
        ax = gca;
        for i = 1:length(I)
            plot(ax,delay(i,:),discAmpl(i,:),'o-','DisplayName',[sel ' = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
            hold on;
            x = delay(i,:);
            if ~any(~isnan(discAmpl(i,:))) %if there are only nans
                discAmpl(i,:) = zeros(1,H);
            end
%             ys = transpose(csaps(x,discAmpl(i,:),0.0001,x));  
%             ys = ys';
            ys = discAmpl(i,:);
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
                    se = fittype(@(a,b,c,d,x) a./(cosh((x-b)/c)).^2+d);                   
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
                    fitPeak(i) = mean(discAmpl(i,Index-2:Index+2)); 
            end;               
        end;
        hold off;
        if strcmp(fitType,'noFit')
            legend('location','best');
        else
            f=get(ax,'Children');
            index = length(f)-((1:length(I))-1).*2;
            legend(f(index),'location','best');
        end
        xlabel(ax,['Delay (' xUnit ')']);
        if plotrelative
            ylabel(ax,'$<Q>/\sqrt{n_{Tg}}$','Interpreter','latex'); 
        else            
            ylabel(ax,'<Q>'); 
        end
        fig = figure(1);
        
    case 'AmplitudeLog'
        figure(1);
        ax = gca;
        for i = 1:length(I)
            semilogy(ax,delay(i,:),discAmpl(i,:)*2^i,'o-','DisplayName',[sel ' = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
            hold on;   
        end
        hold off;
        legend('location','best');
%         if strcmp(fitType,'noFit')
%             legend('location','best');
%         else
%             f=get(ax,'Children');
%             index = length(f)-((1:length(I))-1).*2;
%             legend(f(index),'location','best');
%         end
        xlabel(ax,['Delay (' xUnit ')']); 
        if plotrelative
            ylabel(ax,'$<Q>/\sqrt{n_{Tg}}$','Interpreter','latex'); 
        else            
            ylabel(ax,'<Q>'); 
        end
        fig = figure(1);
        
        
    case 'MeanVar'
        for i = 1:length(I)
            plot(delay(i,:),discMeanVar(i,:),'o-','DisplayName',[sel ' = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
            hold on;
            Index = find(delay(i,:)>= -30 & delay(i,:)<= 30);
            fitPeak(i) = mean(discMeanVar(i,Index)); 
        end;
        hold off;
        legend('location','southeast');
        xlabel(['Delay (' xUnit ')']);
        if plotrelative
            ylabel('Average Variance/n_{Tg}'); 
        else            
            ylabel('Average Variance'); 
        end
     case 'VarQ'
        for i = 1:length(I)
            plot(delay(i,:),varQ(i,:),'o-','DisplayName',[sel ' = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
            hold on;
            Index = find(delay(i,:)>= -30 & delay(i,:)<= 30);
            fitPeak(i) = mean(varQ(i,Index)); 
        end;
        hold off;
        legend('location','southeast');
        if plotrelative
            ylabel('Var_{Q}/n_{Tg}'); 
        else            
            ylabel('Var_{Q}'); 
        end
        xlabel(['Delay (' xUnit ')']);
     case 'VarP'
        for i = 1:length(I)
            plot(delay(i,:),varP(i,:),'o-','DisplayName',[sel ' = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
            hold on;
            Index = find(delay(i,:)>= -30 & delay(i,:)<= 30);
            fitPeak(i) = mean(varP(i,Index)); 
        end;
        hold off;
        legend('location','southeast');
        if plotrelative
            ylabel('Var_{P}/n_{Tg}'); 
        else            
            ylabel('Var_{P}'); 
        end
        xlabel(['Delay (' xUnit ')']);
    case 'MeanVarSigma'
        errorbar(Yr(:,1),sigmas,abs(sigmas-sigmaConf(:,1)), ...
            abs(sigmas-sigmaConf(:,2)),'o-','DisplayName', ...
            'Standard Deviation of Gaussian with 95% confidence intervals');
        xlabel('Ring Radius');
        ylabel('Temporal Width of Minimum Variance');
        title('Width of Minimum Variance vs. Postselected Radius');
        legend('show');
    case 'DiscN'
        figure(1);
        ax = gca;
        for i = 1:length(I)
            plot(ax,delay(i,:),discN(i,:),'o-','DisplayName',[sel ' = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
            hold on;
            x = delay(i,:);
            if ~any(~isnan(discN(i,:))) %if there are only nans
                discN(i,:) = zeros(1,H);
            end
%             ys = transpose(csaps(x,discN(i,:),0.0001,x));  
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
                    se = fittype(@(a,b,c,d,x) a./(cosh((x-b)/c)).^2+d);                   
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
                    fitPeak(i) = mean(discN(i,Index-2:Index+2)); 
            end;               
        end;
        hold off;
        if strcmp(fitType,'noFit')
            legend('location','best');
        else
            f=get(ax,'Children');
            index = length(f)-((1:length(I))-1).*2;
            legend(f(index),'location','best');
        end
        xlabel(ax,['Delay (' xUnit ')']); 
        if plotrelative
            ylabel('n/n_{Tg}'); 
        else            
            ylabel('n'); 
        end
        fig = figure(1);
    case 'G2'
        for i = 1:length(I)
            plot(delay(i,:),g2vals(i,:),'o-','DisplayName',[sel ' = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
            hold on;
        end
        hold off;
        legend('location','northeast');
        xlabel(['Delay (' xUnit ')']);
        ylabel('g^{(2)}(0)');
        ylim([0.5 2.5]);
        if ~varyAPS
            title('G2 vs. Radius of Postselected Fullcircle');
        else
            title('G2 vs. A_c');
        end
    case 'nWithoutPostselection'        
        plot(delay(1,:),nTg,'o-','DisplayName','n_{tg}');
        hold on;
        plot(delay(1,:),nPsFast,'o-','DisplayName','n_{psFast}');
        plot(delay(1,:),nPsSlow,'o-','DisplayName','n_{psSlow}');      
        hold off;
        legend('location','northwest');
        xlabel(['Delay (' xUnit ')']);
        ylabel('mean Photon Number');       
        title('Photon Numbers without Postselection');        
    case 'ThicknessMeanVar'
        surf(delay,Yt,discMeanVar);
        view(-50,20);
        xlabel(['Delay (' xUnit ')']);
        zlabel('Average Variance');
        title('Variance vs. Thickness of Postselected Fullcircle');
    case 'ThicknessMeanVarSigma'
        errorbar(Yt(:,1),sigmas,abs(sigmas-sigmaConf(:,1)), ...
            abs(sigmas-sigmaConf(:,2)),'o','DisplayName', ...
            'Standard Deviation of Gaussian with 95% confidence intervals');
        xlabel('Ring Thickness');
        ylabel('Temporal Width of Minimum Variance');
        title('Width of Minimum Variance vs. Postselected Thickness');
        legend('show');
    case 'ThicknessMeanVarMin'
        [~,iMin] = min(discMeanVar(1,:));
        plot(Yt(:,iMin),discMeanVar(:,iMin));
        xlabel('Ring Thickness');
        ylabel('Minimum of Postselected Variance');
        title('Minimum Variance vs. Thickness of Postselected Fullcircle');
    case 'ThicknessG2'
        surf(delay,Yt,g2vals);
        view(3);
        xlabel(['Delay (' xUnit ')']);
        zlabel('g^{(2)}');
        title('Photon Number vs. Thickness of Postselected Fullcircle');
end
set(fig,'Color','w','Units','centimeters','Position',[1,1,45,30],'PaperPositionMode','auto');
graphicsSettings;
%% Write figure to file and close it
if ~isempty(filename)
    savefig(fig,filename);
    print([filename '.png'],'-dpng');
    close all;
end

%% plot fit results vs ring radius and thickness
if any(fitTau)
    errorbar(Yr(:,1),fitTau,tauError,'o-','Linewidth',2);
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
    savefig([filename '-' typestr '-fitWidthsVsRadius.fig']);
    print([filename '-' typestr '-fitWidthsVsRadius.png'],'-dpng');
    close all;
end

if any(fitPeak)
    plot(Yr(:,1),fitPeak,'o-');
    xlabel('A_{ps} set for postselection');
    ylabel('Peakheight');
    graphicsSettings;
    title(typestr);
    savefig([filename '-' typestr '-PeaksVsRadius.fig']);
    print([filename '-' typestr '-PeaksVsRadius.png'],'-dpng');
    close all;
end

% if any(fitTau)
%     errorbar(Yt(:,1),fitTau,tauError,'o-','Linewidth',2);
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
%     title(typestr);
%     savefig([filename '-' typestr '-fitWidthsVsThickness.fig']);
%     print([filename '-' typestr '-fitWidthsVsThickness.png'],'-dpng');
%     close all;
% end
% 
% if any(fitPeak)
%     plot(Yt(:,1),fitPeak,'o-');
%     xlabel('Ring Thickness set for postselection');
%     ylabel('Peakheight');
%     graphicsSettings;
%     title(typestr);
%     savefig([filename '-' typestr '-PeaksVsThickness.fig']);
%     print([filename '-' typestr '-PeaksVsThickness.png'],'-dpng');
%     close all;
% end
% plot(Yt(:,1),fitWidth,'o-');
% plot(Yt(:,1),fitPeak,'o-');



end

