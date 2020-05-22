function [E0,a,y0,Emax,modeInt, SumInt] = plotDispersion(filenameSIG, filenameBG,varargin)
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
defaultXAperture = 330; 
addParameter(parser,'xAperture',defaultXAperture,@isnumeric);
defaultModeK = [-0.66 0.66]; 
addParameter(parser,'ModeK',defaultModeK,@isnumeric); % k range in which the mode is selected
defaultModeE = [1.611 1.6113]; 
addParameter(parser,'ModeE',defaultModeE,@isnumeric); % Energy range in which the mode is selected
defaultZoomE = [1.605 1.625]; 
addParameter(parser,'ZoomE',defaultZoomE,@isnumeric); % Energy range in which the dispersion is plotted
defaultFit = 'yes'; % if set to 'no', no fit is made. If set to 'useOld', a plot is made with old fit parameters. 
addParameter(parser,'Fit',defaultFit);
defaultE0Old = 1.6108; 
addParameter(parser,'E0Old',defaultE0Old,@isnumeric);
defaultaOld = 4.4973e-7; 
addParameter(parser,'aOld',defaultaOld,@isnumeric);
defaultY0Old = 217; 
addParameter(parser,'y0Old',defaultY0Old,@isnumeric);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[aOld,E0Old,fitoption,minY,modeE,modeK,plotMode,plottype,subtract,xAperture,y0Old,zoomE] = c{:};

%% fourier pixel
% take values from calibration measurement
NA = 0.4; %depends on microscope objective

%% constants
h = 6.62607015e-34;
c0 = 299792458;
e0 = 1.602176634e-19;

%% load data
    cd('raw-data');
    %data = textread(filenameSIG);
    %data = textread(filenameSIG,'','headerlines',1); 
    data = textread(filenameSIG,'','delimiter',',','headerlines',1); 
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
[~,Index] = max(Int(energy<1.63,:));
Emaxs = energy(Index);
minE = Emaxs(minY);
 
%% fit 
if strcmp(fitoption,'yes')  
    % guess value for slope
    aStart = (Emaxs(round(length(Emaxs)*0.75)) - minE)/(y(round(length(y)*0.75))-minY)^2;
    EmaxsFit=Emaxs(Emaxs>=minE & y>= 0.5*minY & y<= 1.5*minY );  %limit the range for the fit 
    yFit=y(Emaxs>=minE & y>= 0.5*minY & y<= 1.5*minY );

    %% make parabolic fit
    par = 'E0 + a * (x-x0)^2';
    [f,gof,~]=fit(yFit',EmaxsFit',par, 'StartPoint', [minE aStart minY]);
    E0 = f.E0;
    a = f.a;
    y0 = f.x0;
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
k = E0*e0*2*pi/(h*c0) * sin(theta);
k = k *1e-6; %k in 1/micrometer

%% get mode integrated intensity and energy of mode maximum
modeRangeE = energy>=min(modeE) & energy <= max(modeE);
modeRangeK = k >=min(modeK) & k <= max(modeK);
modeInt = sum(sum(Int(modeRangeE, modeRangeK)));
modeEnergy = energy(modeRangeE);
[~,index]=max(modeInt);
Emax = modeEnergy(index);
 
%% make 2D surface plot
range = energy>=min(zoomE) & energy <= max(zoomE);
if strcmp(plottype, 'lin')
    pcolor(k, energy(range), Int(range,:));
else
    logInt = log(Int);
    logInt(Int==0) = 0;
    logInt(logInt<0)=0;
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
    plot(k(Efit>=min(zoomE) & Efit <= max(zoomE)),Efit(Efit>=min(zoomE) & Efit <= max(zoomE)),'r-','LineWidth',1.5);
end
if strcmp(fitoption,'useOld')
    Efit = E0Old + aOld * (y-y0Old).^2;
    plot(k(Efit>=min(zoomE) & Efit <= max(zoomE)),Efit(Efit>=min(zoomE) & Efit <= max(zoomE)),'r-','LineWidth',1.5);
end
fontsize = 24;
graphicsSettings;
fontName = 'Times New Roman';
set(gca,'XGrid','on','YGrid','on');
xlabel('k (\mum ^{-1})','FontSize',fontsize,'FontName',fontName);
ylabel('Energy (eV)','FontSize',fontsize,'FontName',fontName);

%write sum
SumInt = sum(sum(Int));
text(0.1, 0.1, ['integr. Int ' num2str(SumInt,'%.0f') ' counts'],'FontSize',fontsize-4,...
    'Units','normalized','Color','w');

cd('..');
print([filenameSIG '-2Dsurfplot-' plottype '-Subtract-' subtract '.png'],'-dpng','-r300');
savefig([filenameSIG '-2Dsurfplot-' plottype '-Subtract-' subtract '.fig']);
clf();

%% make 1D plot 
% %range for integration
% % range = length(x)/8;
% % subInt = Int(:,minX-range:minX+range);
% % intInt = sum(subInt,2);
% 
% % intInt = Int(:,minX);
% % [Intmax,index]=max(intInt);
% % Emax = energy(index);
% intInt = Int(:,minY);
% Intmax = Int(modeY, minY);
% Emax = energy(modeY);
% 
% if strcmp(plottype,'lin')
%     plot(energy, intInt, 'linewidth',2);
%     set(gca, 'Ylim',[0 max(intInt)+10000]);
%     hold on;
% 
% elseif strcmp(plottype,'log')
%     semilogy(energy(intInt>0), intInt(intInt>0), 'linewidth',2);
%     hold on;
%        
% end
% ylabel('Intensity at k = 0 (a.u.)','FontSize',fontsize,'FontName',fontName);
% xlabel('Energy (eV)','FontSize',fontsize,'FontName',fontName);
% set(gca,'LineWidth',3,'XColor',[0 0 0], 'YColor', [0 0 0],'Box','on',...
% 'FontSize',22,'FontName',fontName,...
% 'TickDir','In');
% print([filenameSIG '-k0Plot-' plottype '-Subtract-' subtract '.png'],'-dpng','-r300');
% savefig([filenameSIG '-k0Plot-' plottype '-Subtract-' subtract '.fig']);
% clf();
    
end