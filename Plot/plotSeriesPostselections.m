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
parse(p,varargin{:});
c = struct2cell(p.Results);
[filename,fitType,typestr,varyAPS,xUnit] = c{:};

% Constants
figurepath = 'figures-fig/';

%% Gather data
[delay,Yr,Yt,discAmpl,discMeanVar,discN,g2vals,g2std] = deal([]);
sigmas = zeros(length(listOfParams),1);
sigmaConf = zeros(length(listOfParams),2);
for iParams = 1:length(listOfParams)
    selParams = listOfParams(iParams);
    
    % From tables
    A = seriesRead3ChTable(selParams,'VaryAPS',varyAPS);
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
    discN(iParams,:) = A.discN;
    discN(iParams,:) = discN(iParams,I);
    g2vals(iParams,:) = A.g2;
    g2vals(iParams,:) = g2vals(iParams,I);
    g2std(iParams,:) = A.g2std;
    g2std(iParams,:) = g2std(iParams,I);
    g2std(iParams,:) = A.g2std;
    g2std(iParams,:) = g2std(iParams,I);
    nX1(iParams,:) = A.nX1;
    nX1(iParams,:) = nX1(iParams,I);
    nX2(iParams,:) = A.nX2;
    nX2(iParams,:) = nX2(iParams,I);
    nX3(iParams,:) = A.nX3;
    nX3(iParams,:) = nX3(iParams,I);
    
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
Yr = Yr(I,:); %Radii
Yt = Yt(I,:); %thicknesses
discAmpl = discAmpl(I,:);
discMeanVar = discMeanVar(I,:);
discN = discN(I,:);
g2vals = g2vals(I,:);
g2std = g2std(I,:);
nX1 = nX1(1,:);
nX2 = nX2(1,:);
nX3 = nX3(1,:);
[fitTau,fitPeak,tauError] = deal(zeros(length(I),1));

%% Create figure
fig = figure;
%formatFigA5(fig);
switch typestr
    case 'Amplitude'
        figure(1);
        ax = gca;
        for i = 1:length(I)
            plot(ax,delay(i,:),discAmpl(i,:),'o-','DisplayName',['r = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
            hold on;
            x = delay(i,:);
            if ~any(~isnan(discAmpl(i,:))) %if there are only nans
                discAmpl(i,:) = zeros(1,H);
            end
            %ys = transpose(csaps(x,discAmpl(i,:),0.001,x));  
            ys = discAmpl(i,:);
            switch fitType
                case 'gauss'
                    [res,gof,~] = fit(x',ys,'gauss1'); %f(x) =  a1*exp(-((x-b1)/c1)^2)
                    fitTau(i) = res.c1;
                    level = 2*tcdf(-1,gof.dfe);
                    m = confint(res,level); 
                    tauError(i) = m(end,end) - res.c1;
                    fitPeak(i) = res(res.b1);  
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
                                      
                    initGuess1 = [200,50, 5]; %peak x, gamma, sigma
                    [estimates1, model1] = voigtfit(x, ys, initGuess1, [min(x), max(x)]);
                    res = myvoigt(x, 200, estimates1(2), estimates1(3) );%estimates1(1)
                    res = res';
                    c = abs(res\ys');
                    res = res*c;             
                    plot(ax,x,res,'r','DisplayName','');
            end;               
        end;
        hold off;
        f=get(ax,'Children');
        index = length(f)-((1:length(I))-1).*2;
        legend(f(index),'location','northwest');
        xlabel(['Delay (' xUnit ')']); 
        ylabel('Coherent Amplitude'); 
        if ~varyAPS
            title('Coherent Amplitude vs. Radius of Postselected Fullcircle');
        else
            title('Coherent Amplitude vs. A_c');
        end
        fig = figure(1);
    case 'MeanVar'
        for i = 1:length(I)
            plot(delay(i,:),discMeanVar(i,:),'o-','DisplayName',['r = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
            hold on;
        end;
        hold off;
        legend('location','southeast');
        ylabel('Average Variance');
        xlabel(['Delay (' xUnit ')']);
        if ~varyAPS
            title('Variance vs. Radius of Postselected Fullcircle');
        else
            title('Variance vs. A_c');
        end
    case 'MeanVarSigma'
        errorbar(Yr(:,1),sigmas,abs(sigmas-sigmaConf(:,1)), ...
            abs(sigmas-sigmaConf(:,2)),'o-','DisplayName', ...
            'Standard Deviation of Gaussian with 95% confidence intervals');
        xlabel('Ring Radius');
        ylabel('Temporal Width of Minimum Variance');
        title('Width of Minimum Variance vs. Postselected Radius');
        legend('show');
    case 'DiscN'
        for i = 1:length(I)
            plot(delay(i,:),discN(i,:),'o-','DisplayName',['r = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
            hold on;
        end;
        hold off;
        legend('location','northwest');
        xlabel(['Delay (' xUnit ')']);
        ylabel('Photon Number');
        if ~varyAPS
            title('Photon Number vs. Radius of Postselected Fullcircle');
        else
            title('Photon Number vs. A_c');
        end
    case 'G2'
        for i = 1:length(I)
            plot(delay(i,:),g2vals(i,:),'o-','DisplayName',['r = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
            hold on;
        end;
        hold off;
        legend('location','northeast');
        xlabel(['Delay (' xUnit ')']);
        ylabel('g^{(2)}(0)');
        if ~varyAPS
            title('G2 vs. Radius of Postselected Fullcircle');
        else
            title('G2 vs. A_c');
        end
    case 'nWithoutPostselection'        
        plot(delay(1,:),nX1,'o-','DisplayName','n_1');
        hold on;
        plot(delay(1,:),nX2,'o-','DisplayName','n_2');
        plot(delay(1,:),nX3,'o-','DisplayName','n_3');      
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
set(fig,'Color','w');
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
    xlabel('A_c set for postselection');
    ylabel(['\tau_c of Gaussian (' xUnit ')']);
    graphicsSettings;
    title([filename '-fitWidthsVsRadius']);
    savefig([filename '-fitWidthsVsRadius.fig']);
    print([filename '-fitWidthsVsRadius.png'],'-dpng');
    close all;
end

if any(fitPeak)
    plot(Yr(:,1),fitPeak,'o-');
    xlabel('A_c set for postselection');
    ylabel('Peakheight of Gaussian');
    graphicsSettings;
    title([filename '-fitPeaksVsRadius']);
    savefig([filename '-fitPeaksVsRadius.fig']);
    print([filename '-fitPeaksVsRadius.png'],'-dpng');
    close all;
end

if any(fitTau)
    errorbar(Yt(:,1),fitTau,tauError,'o-','Linewidth',2);
    xlabel('Ring Thickness set for postselection');
    ylabel(['\tau_c of Gaussian (' xUnit ')']);
    graphicsSettings;
    title([filename '-fitWidthsVsThickness']);
    savefig([filename '-fitWidthsVsThickness.fig']);
    print([filename '-fitWidthsVsThickness.png'],'-dpng');
    close all;
end

if any(fitPeak)
    plot(Yt(:,1),fitPeak,'o-');
    xlabel('Ring Thickness set for postselection');
    ylabel('Peakheight of Gaussian');
    graphicsSettings;
    title([filename '-fitPeaksVsThickness']);
    savefig([filename '-fitPeaksVsThickness.fig']);
    print([filename '-fitPeaksVsThickness.png'],'-dpng');
    close all;
end
% plot(Yt(:,1),fitWidth,'o-');
% plot(Yt(:,1),fitPeak,'o-');



end

