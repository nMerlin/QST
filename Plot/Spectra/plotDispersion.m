function [E0,a,y0,Emax,modeInt, SumInt,FWHMFit,FWHMerror,peakPositionFit,peakHeight,Ecav,peakPosition,FWHM,Areapixels] = plotDispersion(filenameSIG, filenameBG,varargin)
%PLOTSTREAK plots the dispersion of polaritons obtained with spectrometer
%and fourier lens and makes a parabolic fit. 
%
%   Input Arguments:
%       filenameSIG: file with the signal data, of 'txt' format.
%       It should be located in a folder 'raw-data'.
%       filenameBG: file with background data.
%       'Subtract','yes': subtract background data
%       minX: supposed position of minimum of parable or middle of k space
%       in pixels
%       xAperture: size of k space in pixels
%
% Output: E0: minimum energy of dispersion. 
% Emax: energy where the maximum of integrated intensity lies. 
% Intmax: maximum of integrated intensity.
% SumInt: Overall Sum of intensity.

%% Validate and parse input arguments
parser = inputParser;
defaultSubtract = 'no'; 
addParameter(parser,'Subtract',defaultSubtract);
defaultPlottype = 'lin'; 
addParameter(parser,'Plottype',defaultPlottype);
defaultPlotMode = 'no'; 
addParameter(parser,'PlotMode',defaultPlotMode); % if yes, the selected mode is indicated in the plot. 
defaultMinY = 200; 
addParameter(parser,'minY',defaultMinY,@isnumeric); % y position of k = 0 (y -->k, x --> energy or wavelength) 
defaultXAperture = 330; %size of full Fourier space in pixels; NA corresponds to half of this. 
addParameter(parser,'xAperture',defaultXAperture,@isnumeric);
defaultModeK = [-0.66 0.66]; 
addParameter(parser,'ModeK',defaultModeK,@isnumeric); % k range in which the mode is selected
defaultModeE = [1.611 1.6113]; 
addParameter(parser,'ModeE',defaultModeE,@isnumeric); % Energy range in which the mode is selected
defaultAdjustEnergy = false; %set whether the energy range of the mode is adjusted to lie around the fitted peak for computing the integrated mode intensity 
addParameter(parser,'AdjustEnergy',defaultAdjustEnergy);
defaultZoomE = [1.59 1.64]; 
addParameter(parser,'ZoomE',defaultZoomE,@isnumeric); % Energy range in which the dispersion is plotted
defaultOnlySmallRange = false;
addParameter(parser,'OnlySmallRange',defaultOnlySmallRange,@islogical); %get FWHM and peakPosition only in the selected energy range
defaultFit = 'yes'; % if set to 'no', no fit is made. If set to 'useOld', a plot is made with old fit parameters. 
addParameter(parser,'Fit',defaultFit);
defaultE0Old = 1.6108; 
addParameter(parser,'E0Old',defaultE0Old,@isnumeric);
defaultaOld = 4.4973e-7; 
addParameter(parser,'aOld',defaultaOld,@isnumeric);
defaultY0Old = 217; 
addParameter(parser,'y0Old',defaultY0Old,@isnumeric);
defaultR = 9.5;  %Rabi splitting in meV
addParameter(parser,'R',defaultR,@isnumeric);
defaultEX = 1.6195;  %Exciton energy in eV
addParameter(parser,'EX',defaultEX,@isnumeric);
defaultPlotOption = true;
addParameter(parser,'PlotOption',defaultPlotOption,@islogical);
defaultNormalize = false; %normalizes the colorscale in the 2D plot to 1.
addParameter(parser,'Normalize',defaultNormalize,@islogical);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[adjustEnergy,aOld,E0Old,EX,fitoption,minY,modeE,modeK,normalize,onlySmallRange,plotMode,plotOption,plottype,R,subtract,xAperture,y0Old,zoomE] = c{:};

%% fourier pixel
% take values from calibration measurement
NA = 0.4; %depends on microscope objective

%% constants
h = 6.62607015e-34;
c0 = 299792458;
e0 = 1.602176634e-19;
nC = 3.7; %refractive index of quantum wells GaAs   

%% load data
    cd('raw-data');
   % data = textread(filenameSIG); % for old measurements with Winspec 
   % data = textread(filenameSIG,'','headerlines',1); 
   data = textread(filenameSIG,'','delimiter',',','headerlines',1);  % for measurements with Lightfield 
    W = data(:,1); % wavelength
    Y = data(:,2); % pixel position
    M = data(:,3); %intensity
    
    %% subtract background
    if strcmp(subtract, 'yes')
        dataBG = textread(filenameBG,'','delimiter',',','headerlines',1);
        IntBG = dataBG(:,3);
        M = M - IntBG;
    else
        %bg = min(M);
        %M = M-bg;
    end
    
    %% set negative values to zero
    M(M<0) = 0;
    
%% sort data
if Y(1) == 0
    Y=Y+1;
end
 n = length(Y)/max(Y);
 w = W(1:n);
 y = 1:max(Y);
 Int = reshape(M, [n max(Y)]);
 Int = flip(Int);
 
%% compute energy in eV from wavelength
 energy = h*c0/e0*1./(w*10^-9);
 energy = fliplr(energy');
 
%% find energy of intensity maxima over the energy axis in the spectral
    % range where we expect the polariton
range = energy>=min(zoomE) & energy <= max(zoomE);
[~,Index] = max(Int(range,:));
E = energy(range);
Emaxs = E(Index);
minE = Emaxs(minY);
 
%% fit 
if strcmp(fitoption,'yes')  
    % guess value for slope
    aStart = (Emaxs(round(1.2*minY)) - minE)/(y(round(1.2*minY))-minY)^2;
    EmaxsFit=Emaxs(Emaxs>=minE-0.05 & y>= 0.69*minY & y<= 1.25*minY );  %limit the range for the fit 
    yFit=y(Emaxs>=minE-0.05 & y>= 0.69*minY & y<= 1.25*minY );


    %% make parabolic fit
    par = 'E0 + a * (x-x0)^2';
    [f,gof,~]=fit(yFit',EmaxsFit',par, 'StartPoint', [minE aStart minY]);
    %[f,gof,~]=fit(yFit',EmaxsFit',par, 'StartPoint', [1.6107 3.61488190561411e-06 210]);
    E0 = f.E0;
    a = f.a;
    y0 = f.x0;   
    %plot for testing:
    %plot(y,Emaxs,'.')
%     plot(yFit',EmaxsFit,'.');
%     plot(y,f(y));
elseif strcmp(fitoption,'useOld')
    E0 = E0Old;
    y0 = y0Old;
    a = aOld;    
else
    E0 = minE;
    y0 = minY;
    a = 0;
end
 
%% compute k from y
theta_max = asin(NA);
theta = (y - y0)*2*theta_max/xAperture;
k = E0*e0*2*pi/(h*c0) * sin(theta); %See Kasprak 2006:Bose?Einstein condensation ofexciton polaritons; 
%E0 is used as an approximation for E.
k = k *1e-6; %k in 1/micrometer

%% get mode integrated intensity and energy of mode maximum
modeRangeE = energy>=min(modeE) & energy <= max(modeE);
modeRangeK = k >=min(modeK) & k <= max(modeK);
modeIntPreliminary = sum(sum(Int(modeRangeE, modeRangeK)));
modeEnergy = energy(modeRangeE);
[~,index]=max(modeIntPreliminary);
Emax = modeEnergy(index);

%get in k direction integrated energy, maximum and peak position and FWHM for later use 
subInt = Int(range, modeRangeK);
intInt = sum(subInt,2);
intInt = intInt - mean(intInt(1:20)); %subtract background 
zoomEnergy = energy(range);
[Max,I] = max(intInt);
peakPosition = zoomEnergy(I);
%get integrated intensity around mode; if wished, adjust the energy range
%of the mode according to the found peak position 
if adjustEnergy
    modeRangeE = energy>= peakPosition - 4e-4 & energy <= peakPosition + 4e-4;
    modeEnergy = energy(modeRangeE);
end
modeInt = sum(sum(Int(modeRangeE, modeRangeK))); % integrated intensity in the selected mode rectangle. 

%get FWHM and peakPosition only in the selected energy range, if chosen
if onlySmallRange
    intInt2 = Int(:, modeRangeK);
    intInt2 = sum(intInt2,2);
    intInt2 = intInt2 - mean(intInt2(1:20));
    try
        FWHM = fwhm(modeEnergy,intInt2(modeRangeE));
    catch
        FWHM = fwhm(zoomEnergy, intInt)* 1000; %in meV;
    end
    [~,I2] = max(intInt2(modeRangeE));
    peakPosition = modeEnergy(I2);
else
    FWHM = fwhm(zoomEnergy, intInt)* 1000; %in meV;
end

A = Int(modeRangeE, modeRangeK);
Areapixels = length(A(:)); % number of pixels in the selected mode rectangle. 
 
%% make 2D surface plot
logInt = log(Int);
logInt(Int==0) = 0;
logInt(logInt<0)=0;
if strcmp(plottype, 'lin')
    if normalize
        IntNorm = Int/max(max(Int));
        pcolor(k, energy(range), IntNorm(range,:));
    else
        pcolor(k, energy(range), Int(range,:));
    end
else
    pcolor(k, energy(range), logInt(range,:));
end
colorbar;
shading flat;
axis tight;
hold on;
if strcmp(plotMode,'yes')
   pcolor(k(modeRangeK),modeEnergy,logInt(modeRangeE, modeRangeK));
end
if strcmp(fitoption,'yes')
    Efit = f(y);
end
if strcmp(fitoption,'useOld')
    Efit = E0Old + aOld * (y-y0Old).^2;
end
if strcmp(fitoption,'yes') || strcmp(fitoption,'useOld')
    [detuningMeV, detuning2g0,Ecav,~] = computeDetuning(E0,R,EX(1));
    EX = EX*ones(length(k(Efit>=min(zoomE) & Efit <= max(zoomE))),1);   
    %Ecav = Ecav*ones(length(k(Efit>=min(zoomE) & Efit <= max(zoomE))),1);
    mcav = Ecav*e0 *nC^2 /c0^2;
    EcavParabola = Ecav + (h/(2*pi))^2 * (k(Efit>=min(zoomE) & Efit <= max(zoomE))*1e6).^2 ./(2*mcav) /e0;
    plot(k(Efit>=min(zoomE) & Efit <= max(zoomE)),Efit(Efit>=min(zoomE) & Efit <= max(zoomE)),...
        'r-','LineWidth',1.5,'DisplayName','E_{LP}');
    plot(k(Efit>=min(zoomE) & Efit <= max(zoomE)),EcavParabola,'y--','LineWidth',1.5,'DisplayName','E_{cav}');
    plot(k(Efit>=min(zoomE) & Efit <= max(zoomE)),EX,'w-','LineWidth',1.5,'DisplayName','E_{exc}');
end
ax = gca;
graphicsSettings;
graphs=get(ax,'Children');
legend(graphs(1:end-1),'Location','best','TextColor','w'); 
legend('boxoff');
fontsize = 24;
fontName = 'Arial';
set(gca,'XGrid','on','YGrid','on');
xlim([-3,3]);
ylim(zoomE);
xlabel('k (\mum ^{-1})','FontSize',fontsize,'FontName',fontName);
ylabel('Energy (eV)','FontSize',fontsize,'FontName',fontName);

%write sum
SumInt = sum(sum(Int));
% text(0.1, 0.1, ['integr. Int ' num2str(SumInt,'%.0f') ' counts'],'FontSize',fontsize-4,...
%     'Units','normalized','Color','w');

cd('..');
if plotOption
figure = gcf;
figure.InvertHardcopy = 'off'; 
figure.Color = 'w';
print([filenameSIG '-2Dsurfplot-' plottype '-Subtract-' subtract '-PlotMode-' ...
    plotMode '-adjustModeEnergy-' num2str(adjustEnergy) '-modeE-' num2str(min(modeE)) '-' num2str(max(modeE)) ...
     '-modeK-' num2str(min(modeK)) '-' num2str(max(modeK)) '-normalize-' num2str(normalize) '.png'],'-dpng','-r300');
savefig([filenameSIG '-2Dsurfplot-' plottype '-Subtract-' subtract '-PlotMode-' ...
    plotMode '-adjustModeEnergy-' num2str(adjustEnergy) '-modeE-' num2str(min(modeE)) '-' num2str(max(modeE)) ...
     '-modeK-' num2str(min(modeK)) '-' num2str(max(modeK)) '-normalize-' num2str(normalize) '.fig']);
end
clf();

%% make 1D plot 

plot(zoomEnergy, intInt, 'linewidth',2);
%set(gca, 'Ylim',[0 max(intInt)+10000]);
hold on;   
% make Gauss Fit 
gaussCustom = 'a1*exp(-((x-b1)/c1)^2)';
wfit = zoomEnergy';
Intfit = intInt;
[f,gof,~] = fit(wfit,Intfit,gaussCustom, 'StartPoint', [Max, peakPosition, 0.001] ); %f(x) =  a1*exp(-((x-b1)/c1)^2)
peakPositionFit = f.b1;
peakHeight = f.a1;
FWHMFit = 2*f.c1*sqrt(log(2));
level = 2*tcdf(-1,gof.dfe);
m = confint(f,level); 
std = m(end,end) - f.c1;
FWHMerror = std * 2*sqrt(log(2));
FWHMFit = FWHMFit * 1000; %in meV
FWHMerror = FWHMerror * 1000;
energyPlot = min(zoomEnergy):0.00001:max(zoomEnergy);
plot(energyPlot,f(energyPlot),'r','LineWidth',2,'DisplayName','Fit');
set(gca,'DefaultTextInterpreter','latex');
ylim([0 Max*1.1]);
text(peakPositionFit-0.05, Max/2,...
    ['FWHM = ' num2str(FWHMFit,'%.3f') ' $\pm$ ' num2str(FWHMerror,'%.4f') ' meV' ],'FontSize',fontsize-4);
ylabel('integrated Intensity around k = 0 (a.u.)','FontSize',fontsize,'FontName',fontName);
xlabel('Energy (eV)','FontSize',fontsize,'FontName',fontName);
graphicsSettings;
if plotOption
print([filenameSIG '-modeE-' num2str(min(modeE)) '-' num2str(max(modeE)) ...
     '-modeK-' num2str(min(modeK)) '-' num2str(max(modeK)) '-cut.png'],'-dpng','-r300');
savefig([filenameSIG '-modeE-' num2str(min(modeE)) '-' num2str(max(modeE)) ...
     '-modeK-' num2str(min(modeK)) '-' num2str(max(modeK)) '-cut.fig']);
end
clf();
    
end