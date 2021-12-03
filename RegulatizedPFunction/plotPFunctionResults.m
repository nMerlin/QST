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
defaultRemoveBaseline = false; 
addParameter(p,'RemoveBaseline',defaultRemoveBaseline,@islogical);
defaultZeroDelay = 0;
addParameter(p,'ZeroDelay',defaultZeroDelay,@isnumeric); % in mm
defaultPlotrelative = false; %plots the data relative to the target photon number
addParameter(p,'Plotrelative',defaultPlotrelative,@islogical);
defaultLogplot = 'false'; %makes yscale of plot logarithmic in case 'true'; in case 'shifted', the lines will be stacked to see them better
addParameter(p,'Logplot',defaultLogplot,@isstr);
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
[fitType,logplot,maxQuad,maxX,norm,phiStep,plotrelative,range,removeBaseline,removeMod,res,rvalue,varyAPS,XStep,xUnit,zeroDelay] = c{:};

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
[delay,Yr,Yt,PhotonNrs,meanPh,meanR,varPh,circVa1,circVa2,circVa1Err,circVa2Err,varR,Rerr,phErr,varRerr,varPhErr,PhotonNrErrs] = deal([]);
for iParams = 1:length(listOfParams)
    selParams = listOfParams(iParams);
    selStr = selParamsToStr(selParams);
    foldername = ['Pfunctionplots-',selStr,'-remMod-',...
            num2str(removeMod),'-range-',num2str(range),'-varyAPS-',num2str(varyAPS),'-R-' num2str(rvalue) ...
        '-maxQuad-' num2str(maxQuad) '-Resolution-' num2str(res) '-maxX-' num2str(maxX) '-Xstep-' num2str(XStep) '-phiStep-' num2str(phiStep)];
    load([foldername '\Pfunctionresults.mat'],'Delay','DelayMm','Pmax','sigmaPmax','meanPhase','meanPhaseBinned','varPhase',...
    'varPhaseBinned','g1','sigNeg','meanAmp','meanAmpBinned','circVar1','circVar2','varAmp','varAmpBinned','PhotonNr','PhotonNrBinned',...
    'phaseErr','ampErr','varPhaseErr','varAmpErr','PhotonNrErr','circVar1Err','circVar2Err');
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
    circVa1(iParams,:) = circVar1;
    circVa2(iParams,:) = circVar2;
    varR(iParams,:) = varAmp;
    Rerr(iParams,:) = ampErr;
    phErr(iParams,:) = phaseErr;
    varRerr(iParams,:) = varAmpErr;
    varPhErr(iParams,:) = varPhaseErr;
    PhotonNrErrs(iParams,:) = PhotonNrErr;  
    circVa1Err(iParams,:) = circVar1Err;
    circVa2Err(iParams,:) = circVar2Err;
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
meanR = meanR(I,:);varPh = varPh(I,:);varR = varR(I,:);circVa1 = circVa1(I,:); circVa2 = circVa2(I,:);
circVa1Err = circVa1Err(I,:); circVa2Err = circVa2Err(I,:);
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
%typestrVector = {'R','discN','varR','meanPh','varPh','RLog'};
%typestrVector = {'R'};
typestrVector = {'circVa1'};
for typeI = 1:length(typestrVector)
    plotStuff(cell2mat(typestrVector(typeI)))
end

    function [] = plotStuff(typestr)
%% Create figure
fig = figure;
filename = [figurepath,typestr,'-',fitType,'-remMod-',...
            num2str(removeMod),'-range-',num2str(min(range)),...
            '-',num2str(max(range)),'-varyAPS-',num2str(varyAPS),...
            '-plotrelative-' num2str(plotrelative),...
            '-logplot-',num2str(logplot),'-remBaseline-' num2str(removeBaseline), '.fig'];
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
        case 'circVa1'
            ys = circVa1(i,:); 
            yErr = circVa1Err(i,:); 
             if strcmp(logplot,'true') || strcmp(logplot,'shifted')
                 ylabel(ax,'1 - Circ. Var \phi');
             else                 
                ylabel(ax,'Circ. Var \phi');
             end
        case 'circVa2'
            ys = circVa2(i,:); 
            yErr = circVa2Err(i,:); 
            ylabel(ax,'Circular variance \phi 2');
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

    if strcmp(logplot,'shifted')
        % for circ. varPhi
        shadedErrorBar(delay(i,:),exp(i)*(1-ys),yErr,'lineProps',{'o-','DisplayName',[sel ' = ' num2str(Yr(i,1),2) ', t = ' num2str(Yt(i,1),2)]});
    elseif strcmp(logplot,'true')
        shadedErrorBar(delay(i,:),1-ys,yErr,'lineProps',{'o-','DisplayName',[sel ' = ' num2str(Yr(i,1),2) ', t = ' num2str(Yt(i,1),2)]});
    else
        shadedErrorBar(delay(i,:),ys,yErr,'lineProps',{'o-','DisplayName',[sel ' = ' num2str(Yr(i,1),2) ', t = ' num2str(Yt(i,1),2)]});
    end
    %plot(ax,delay(i,:),ys,'o-','DisplayName',[sel ' = ' num2str(Yr(i,1)) ', t = ' num2str(Yt(i,1))]);
    hold on;
    x = delay(i,:);
    if ~any(~isnan(ys)) %if there are only nans
        ys = zeros(1,H);
    end
    
    %remove nans
    x = x(~isnan(ys));
    ys = ys(~isnan(ys));
    
    if removeBaseline
        z = smooth(ys',400);
        z = z';
        ys = ys - z + mean(z) ;
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
        case 'gaussSat1'
            g = fittype(@(a,b,c,x) a*exp(-pi/2*((x-b)/c).^2)+1);                   
            b0 = zeroDelay;
            c0 = 0.5*max(x);
            d0 = mean(ys(end-5:end)); % saturation value
            a0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay))) - d0;
            [result,gof,~] = fit(x',ys',g,'StartPoint',[a0 b0 c0]); 
            [se]= getStandardErrorsFromFit(result,gof,'method1'); 
            pa1(i) = result.a;
            pa2(i) = result.b;
            pa3(i) = result.c;
            pa1Err(i) = se(1);
            pa2Err(i) = se(2);
            pa3Err(i) = se(3);
            fitTau(i) = result.c;            
            tauErr(i) = se(3);
            fitPeak(i) = result.a + 1; 
            fitPeakErr(i) = sqrt(se(1)^2);           
            plot(ax,min(delay(i,:)):1:max(delay(i,:)),result(min(delay(i,:)):1:max(delay(i,:))),'r','DisplayName',''); 
        case 'envelopeFromPeaks'
            b0 = zeroDelay;
            c0 = 0.5*max(x);
            d0 = mean(ys(end-5:end)); % saturation value
            a0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay))) - d0;
            if a0 < 0
                [peaks,Index] = findpeaks(-ys,'MinPeakProminence',0.02);
                 peaks = -peaks;
            else
                [peaks,Index] = findpeaks(ys,'MinPeakProminence',0.02);
            end
            pDelays = x(Index);
            xFit = pDelays';
            yFit = peaks';
            g = fittype(@(a,b,c,d,x) a*exp(-pi/2*((x-b)/c).^2)+d);                              
            [result,gof,~] = fit(xFit,yFit,g,'StartPoint',[a0 b0 c0 d0]); 
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
         case 'envelope'
            b0 = zeroDelay;
            c0 = 0.5*max(x);
            d0 = mean(ys(end-5:end)); % saturation value
            a0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay))) - d0;
            [yupper,ylower] = envelope(ys,10,'peak');
            if a0 < 0
                yFit = ylower;
            else
                yFit = yupper;
            end
            % for envelope, only small range of delay is important
            xFit = x(delay(i,:)>= (-600+zeroDelay) & delay(i,:)<= (600+zeroDelay));
            yFit = yFit(delay(i,:)>= (-600+zeroDelay) & delay(i,:)<= (600+zeroDelay));            
            g = fittype(@(a,b,c,x) a*exp(-pi/2*((x-b)/c).^2)+d0);                              
            [result,gof,~] = fit(xFit',yFit',g,'StartPoint',[a0 b0 c0]); 
            [se]= getStandardErrorsFromFit(result,gof,'method1'); 
            pa1(i) = result.a;
            pa2(i) = result.b;
            pa3(i) = result.c;
            %pa4(i) = result.d;
            pa1Err(i) = se(1);
            pa2Err(i) = se(2);
            pa3Err(i) = se(3);
            %pa4Err(i) = se(4);
            fitTau(i) = result.c;            
            tauErr(i) = se(3);
%             fitPeak(i) = result.a + result.d; 
%             fitPeakErr(i) = sqrt(se(1)^2 + se(4)^2); 
            fitPeak(i) = result.a + d0; 
            fitPeakErr(i) = sqrt(se(1)^2);   
            plot(ax,min(delay(i,:)):1:max(delay(i,:)),result(min(delay(i,:)):1:max(delay(i,:))),'r','DisplayName','');
        case 'envelopeExp2'
            b0 = zeroDelay;
            c0 = 0.5*max(x);
            d0 = mean(ys(end-5:end)); % saturation value
            a0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay))) - d0;
            [yupper,ylower] = envelope(ys,10,'peak');
            if a0 < 0
                yFit = ylower;
            else
                yFit = yupper;
            end
            g = fittype(@(a,b,c,x) a*exp(-((x-b)/c)) + d0 );
            % for envelope, only small range of delay is important
            xFit = x;
            xFit = x(delay(i,:)>= (-600+zeroDelay) & delay(i,:)<= (600+zeroDelay));
            yFit = yFit(delay(i,:)>= (-600+zeroDelay) & delay(i,:)<= (600+zeroDelay));
            [result,gof,~] = fit(abs(xFit)',yFit',g,'StartPoint',[a0 b0 c0]); 
            [se]= getStandardErrorsFromFit(result,gof,'method1'); 
            pa1(i) = result.a;
            pa2(i) = result.b;
            pa3(i) = result.c;
            %pa4(i) = result.d;
            pa1Err(i) = se(1);
            pa2Err(i) = se(2);
            pa3Err(i) = se(3);
            %pa4Err(i) = se(4);
            fitTau(i) = result.c;            
            tauErr(i) = se(3);
            %fitPeak(i) = result.a + result.d; 
            fitPeak(i) = result.a + d0; 
            %fitPeakErr(i) = sqrt(se(1)^2 + se(4)^2);     
            fitPeakErr(i) = sqrt(se(1)^2);   
            plot(ax,min(delay(i,:)):1:max(delay(i,:)),result(abs(min(delay(i,:)):1:max(delay(i,:)))),'r','DisplayName','');                    
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
            if strcmp(logplot,'true')
                plot(ax,min(delay(i,:)):1:max(delay(i,:)),1-result(abs(min(delay(i,:)):1:max(delay(i,:)))),'r','DisplayName','');
            elseif strcmp(logplot,'shifted')
                plot(ax,min(delay(i,:)):1:max(delay(i,:)),exp(i)*(1-result(abs(min(delay(i,:)):1:max(delay(i,:))))),'r','DisplayName','');
            else
                plot(ax,min(delay(i,:)):1:max(delay(i,:)),result(abs(min(delay(i,:)):1:max(delay(i,:)))),'r','DisplayName','');
            end
        case 'exp2sat1'
            g = fittype(@(a,b,c,x) a*exp(-((x-b)/c)) + 1 );                   
            b0 = zeroDelay;
            c0 = 0.5*max(x);
            d0 = mean(ys(end-5:end)); % saturation value
            a0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay))) - d0;
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
            fitPeak(i) = result.a + 1; 
            fitPeakErr(i) = sqrt(se(1)^2); 
            if strcmp(logplot,'true')
                plot(ax,min(delay(i,:)):1:max(delay(i,:)),1-result(abs(min(delay(i,:)):1:max(delay(i,:)))),'r','DisplayName','');
                residuals = (1-ys') - (1-result(abs(x)));
            elseif strcmp(logplot,'shifted')
                plot(ax,min(delay(i,:)):1:max(delay(i,:)),exp(i)*(1-result(abs(min(delay(i,:)):1:max(delay(i,:))))),'r','DisplayName','');
                residuals = (1-ys') - (1-result(abs(x)));
            else
                plot(ax,min(delay(i,:)):1:max(delay(i,:)),result(abs(min(delay(i,:)):1:max(delay(i,:)))),'r','DisplayName','');
                residuals = ys' - result(abs(x));
            end
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
        case 'power-law'
            g = fittype(@(a,alpha,x) a*abs(x).^(-alpha));                   
            b0 = zeroDelay;
            alpha0 = 0.5;
            a0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay)));
            [result,gof,~] = fit( x(abs(x) >= (100+zeroDelay))',...
                 ys(abs(x)>= (100+zeroDelay))',g,'StartPoint',[a0 alpha0]); 
                        [se]= getStandardErrorsFromFit(result,gof,'method1'); 
            pa1(i) = result.a;
            pa2(i) = result.alpha;
            pa1Err(i) = se(1);
            pa2Err(i) = se(2);  
            fitTau(i) = result.alpha;
            tauErr(i) = se(1);
            fitPeak(i) = result.a; 
            fitPeakErr(i) = se(1);  
            % For the power law, there is no coherence time, because an
            % integral over it doesnt converge. 
            if strcmp(logplot,'true')
                plot(ax,[min(delay(i,:)):1:-100+zeroDelay 100+zeroDelay:1:max(delay(i,:))],...
                (1-result(abs([min(delay(i,:)):1:-100+zeroDelay 100+zeroDelay:1:max(delay(i,:))]))),'r','DisplayName','');  
            elseif strcmp(logplot,'shifted')
                plot(ax,[min(delay(i,:)):1:-100+zeroDelay 100+zeroDelay:1:max(delay(i,:))],...
                exp(i)*(1-result(abs([min(delay(i,:)):1:-100+zeroDelay 100+zeroDelay:1:max(delay(i,:))]))),'r','DisplayName','');   
            else
               plot(ax,[min(delay(i,:)):1:-100+zeroDelay 100+zeroDelay:1:max(delay(i,:))],...
                result(abs([min(delay(i,:)):1:-100+zeroDelay 100+zeroDelay:1:max(delay(i,:))])),'r','DisplayName','');   
            end
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
             if strcmp(logplot,'true')
                plot(ax,min(delay(i,:)):1:max(delay(i,:)),result(abs(min(delay(i,:)):1:max(delay(i,:)))),'r','DisplayName','');
            elseif strcmp(logplot,'shifted')
                plot(ax,min(delay(i,:)):1:max(delay(i,:)),exp(i)*(result(abs(min(delay(i,:)):1:max(delay(i,:))))),'r','DisplayName','');
            else
                plot(ax,min(delay(i,:)):1:max(delay(i,:)),result(abs(min(delay(i,:)):1:max(delay(i,:)))),'r','DisplayName','');
             end   
         case 'stretched-exp-sat1'
            g = fittype(@(A,B,beta,x) A*exp(-B*abs(x).^beta) + 1);                   
            B0 = 1/max(x);
            beta0 = 1;
            A0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay)))-1;
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
             if strcmp(logplot,'true')
                plot(ax,min(delay(i,:)):1:max(delay(i,:)),1-result(abs(min(delay(i,:)):1:max(delay(i,:)))),'r','DisplayName','');
            elseif strcmp(logplot,'shifted')
                plot(ax,min(delay(i,:)):1:max(delay(i,:)),exp(i)*(1-result(abs(min(delay(i,:)):1:max(delay(i,:))))),'r','DisplayName','');
            else
                plot(ax,min(delay(i,:)):1:max(delay(i,:)),result(abs(min(delay(i,:)):1:max(delay(i,:)))),'r','DisplayName','');
             end   
         case 'power-law-1-'
            g = fittype(@(a,alpha,x) a*abs(x).^(-alpha));                   
            b0 = zeroDelay;
            alpha0 = 0.5;
            a0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay)));
            [result,gof,~] = fit( x(abs(x)>= (100+zeroDelay))',...
                 1-ys(abs(x)>= (100+zeroDelay))',g,'StartPoint',[a0 alpha0],'Upper',[1 1]); 
                        [se]= getStandardErrorsFromFit(result,gof,'method1'); 
            pa1(i) = result.a;
            pa2(i) = result.alpha;
            pa1Err(i) = se(1);
            pa2Err(i) = se(2);  
            fitTau(i) = result.alpha;
            tauErr(i) = se(1);
            fitPeak(i) = result.a; 
            fitPeakErr(i) = se(1);  
            % For the power law, there is no coherence time, because an
            % integral over it doesnt converge. 
            if strcmp(logplot,'true')
                plot(ax,min(delay(i,:)):1:max(delay(i,:)),...
                (result(abs(min(delay(i,:)):1:max(delay(i,:))))),'r','DisplayName','');  
                residuals = (1-ys') - result(abs(x));
            elseif strcmp(logplot,'shifted')
                plot(ax,min(delay(i,:)):1:max(delay(i,:)),...
                exp(i)*(result(abs(min(delay(i,:)):1:max(delay(i,:))))),'r','DisplayName','');  
                residuals = (1-ys') - result(abs(x));
            else
               plot(ax,min(delay(i,:)):1:max(delay(i,:)),...
                1-result(abs(min(delay(i,:)):1:max(delay(i,:)))),'r','DisplayName','');   
                residuals = (ys') - (1-result(abs(x)));
            end
       case 'power-law-1-shifted'
            g = fittype(@(a,b,alpha,x) real(a*(abs(x)-b).^(-alpha)));                   
            b0 = zeroDelay;
            alpha0 = 0.5;
            a0 = mean(ys(delay(i,:)>= (-30+zeroDelay) & delay(i,:)<= (30+zeroDelay)));
            [result,gof,~] = fit( x',...
                 1-ys',g,'StartPoint',[a0 b0 alpha0],'Upper',[1 max(x) 1],'Lower',[0 0 0]); 
                        [se]= getStandardErrorsFromFit(result,gof,'method1'); 
            pa1(i) = result.a;
            pa2(i) = result.alpha;
            pa3(i) = result.b;
            br = result.b;
            pa1Err(i) = se(1);
            pa2Err(i) = se(2);  
            pa3Err(i) = se(3); 
            fitTau(i) = result.alpha;
            tauErr(i) = se(1);
            fitPeak(i) = result.a; 
            fitPeakErr(i) = se(1);  
            % For the power law, there is no coherence time, because an
            % integral over it doesnt converge. 
            if strcmp(logplot,'true')
                plot(ax,[min(delay(i,:)):1:-br+zeroDelay br+zeroDelay:1:max(delay(i,:))],...
                (result(abs([min(delay(i,:)):1:-br+zeroDelay br+zeroDelay:1:max(delay(i,:))]))),'r','DisplayName','');
                residuals = (1-ys') - result(abs(x));
            elseif strcmp(logplot,'shifted')
                plot(ax,[min(delay(i,:)):1:-br+zeroDelay br+zeroDelay:1:max(delay(i,:))],...
                exp(i)*(result(abs([min(delay(i,:)):1:-br+zeroDelay br+zeroDelay:1:max(delay(i,:))]))),'r','DisplayName','');   
                residuals = (1-ys') - result(abs(x));
            else
               plot(ax,[min(delay(i,:)):1:-br+zeroDelay br+zeroDelay:1:max(delay(i,:))],...
                1-result(abs([min(delay(i,:)):1:-br+zeroDelay br+zeroDelay:1:max(delay(i,:))])),'r','DisplayName','');  
                residuals = (ys') - (1-result(abs(x)));
            end  
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
    if strcmp(logplot,'true') || strcmp(logplot,'shifted')
        l = legend('location','bestoutside');
    else
        l = legend('location','southeast');
    end
    l.FontSize = 30;
else
    f=get(ax,'Children');
    index = length(f)-((1:length(I))-1).*2;
     if strcmp(logplot,'true') || strcmp(logplot,'shifted')
        l = legend(f(index),'location','bestoutside');
    else        
        l = legend(f(index),'location','southeast');
    end
    l.FontSize = 30;
end
xlabel(ax,['Time delay \tau (' xUnit ')']);
 if strcmp(logplot,'true') || strcmp(logplot,'shifted')
    ax.YScale = 'log';
end
fig = figure(1);              
set(fig,'Color','w','Units','centimeters','Position',[1,1,45,30],'PaperPositionMode','auto');
graphicsSettings;
set(ax,'FontSize',42);
%% Write figure to file and close it
if ~isempty(filename)
    savefig(fig,filename);
    print([filename '.png'],'-dpng');
    close all;
end

%% plot fit results vs ring radius and thickness
if any(fitTau)
    YrPlot = Yr(:,1);
    Yt = Yt(:,1);
    [YrPlot,Ir]= sort(YrPlot);
    Yt = Yt(Ir);
    fitTau = real(fitTau(Ir));
    tauErr = tauErr(Ir);
    errorbar(YrPlot,fitTau,tauErr,'o-','Linewidth',2);
    xlabel('r_{ps} set for postselection');
    switch fitType
        case 'envelopeExp2'
            ylabel(['\tau_c of exp. function (' xUnit ')']);
        case 'gauss'
            ylabel(['\tau_c of Gaussian (' xUnit ')']);
        case 'gaussSat1'
            ylabel(['\tau_c of Gaussian (' xUnit ')']);
        case 'envelope'
            ylabel(['\tau_c of Gaussian (' xUnit ')']);
        case 'voigt'
            ylabel(['FWHM of Voigt Profile (' xUnit ')']);
        case 'exponential'
            ylabel(['\tau_c of exp. function (' xUnit ')']);
        case 'exp2'
            ylabel(['\tau_c of exp. function (' xUnit ')']);
        case 'exp3'
            ylabel(['\tau_c of exp. function (' xUnit ')']);
        case 'exp2sat1'
            ylabel(['\tau_c of exp. function (' xUnit ')']);
        case 'sech'
            ylabel(['\tau_c of sech function (' xUnit ')']);
         case 'lorentz'
            ylabel(['FWHM of Lorentz Profile (' xUnit ')']);
         case {'power-law','power-law-1-'}
            ylabel('\alpha of power-law Profile');
        case 'stretched-exp'
            ylabel('\tau_c of stretched-exp Profile');
    end
    graphicsSettings;
    ax = gca;
    set(ax,'FontSize',25);
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
        case {'power-law','power-law-1-'}
            ylabel('A of power-law Profile');
        otherwise
            ylabel('Peakheight');
    end
    graphicsSettings;
    ax = gca;
    set(ax,'FontSize',25);
    title(typestr);
    savefig([filename '-PeaksVsRadius.fig']);
    print([filename '-PeaksVsRadius.png'],'-dpng');
    close all;
end

%% plot fit residuals
if any(residuals)
    plot(x,residuals);
    hold on;
    plot(x,zeros(length(x)),'k--');
    xlabel(['Time delay \tau (' xUnit ')']);
    ylabel('Residuals for highest r_{ps}');
    ylim([-0.5 0.5]);
    graphicsSettings;
    ax = gca;
    set(ax,'FontSize',24);
    title(typestr);
    savefig([filename '-Residuals.fig']);
    print([filename '-Residuals.png'],'-dpng');
    close all;
end

%% save fitresults
save([filename '-fitresults.mat'],'YrPlot','Yt','fitTau','tauErr','fitPeak','fitPeakErr',...
    'sse','rsquare','adjrsquare','rmse','pa1','pa2','pa3','pa4','pa1Err','pa2Err','pa3Err','pa4Err');
T = table(YrPlot,Yt,fitTau,tauErr,fitPeak,fitPeakErr,sse,rsquare,adjrsquare,rmse,pa1,pa2,pa3,pa4,pa1Err,pa2Err,pa3Err,pa4Err);
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

