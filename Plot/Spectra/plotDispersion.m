function [E0] = plotDispersion(filenameSIG, filenameBG,varargin)
%PLOTSTREAK plots the dispersion of polaritons obtained with spectrometer
%and fourier lens and makes a parabolic fit. 
%
%   Input Arguments:
%       filenameSIG: file with the signal data, of 'txt' format.
%       It should be located in a folder 'raw-data'.
%       filenameBG: file with background data.
%       'Subtract','yes': subtract background data
%
% Output: E0: minimum energy of dispersion. 

%% Validate and parse input arguments
parser = inputParser;
defaultSubtract = 'no'; 
addParameter(parser,'Subtract',defaultSubtract);
defaultPlottype = 'lin'; 
addParameter(parser,'Plottype',defaultPlottype);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[plottype,subtract] = c{:};

%% fourier pixel
% take values from calibration measurement
NA = 0.4; %depends on microscope objective
minX = 202;
xAperture = 330;

%% constants
h = 6.62607015e-34;
c0 = 299792458;
e0 = 1.602176634e-19;

%% load data
    cd('raw-data');
    data = textread(filenameSIG);
    W = data(:,1); % wavelength
    X = data(:,2); % pixel position
    M = data(:,3); %intensity
    
    %% subtract background
    if strcmp(subtract, 'yes')
        dataBG = textread(filenameBG);
        IntBG = dataBG(:,3);
        M = M - IntBG;
    else
        %bg = min(M);
        %M = M-bg;
    end
    
%% sort data
 n = length(X)/max(X);
 w = W(1:n);
 x = 1:max(X);
 Int = reshape(M, [n max(X)]);
 Int = flip(Int);
 
%% compute energy in eV from wavelength
 energy = h*c0/e0*1./(w*10^-9);
 energy = fliplr(energy');

% find energy of intensity maxima over the energy axis in the spectral
% range where we expect the polariton
[~,Index] = max(Int(energy<1.63,:));
Emaxs = energy(Index);
minE = Emaxs(minX);
% guess value for slope
aStart = (Emaxs(round(length(Emaxs)*0.75)) - minE)/(x(round(length(x)*0.75))-minX)^2;
Emaxs=Emaxs(Emaxs>=minE);
xFit=x(Emaxs>=minE);
    
%% make parabolic fit
par = 'E0 + a * (x-x0)^2';
[f,gof,~]=fit(xFit',Emaxs',par, 'StartPoint', [minE aStart minX]);
E0 = f.E0;
a = f.a;
x0 = f.x0;
 
%% compute k from x
theta_max = asin(NA);
theta = (x - minX)*2*theta_max/xAperture;
k = E0*e0*2*pi/(h*c0) * sin(theta);
k = k *1e-6; %k in 1/micrometer
 
%% make 2D surface plot
if strcmp(plottype, 'lin')
    surf(k, energy, Int);
else
    surf(k, energy, log(Int));
end
colorbar;
view(180,-90);
shading flat;
axis tight;
hold on;
plot(k,f(x),'r-','LineWidth',1.5);

fontsize = 24;
graphicsSettings;
fontName = 'Times New Roman';
set(gca,'XGrid','on','YGrid','on');
xlabel('k (\mum ^{-1})','FontSize',fontsize,'FontName',fontName);
ylabel('Energy (eV)','FontSize',fontsize,'FontName',fontName);

cd('..');
print([filenameSIG '-2Dsurfplot.png'],'-dpng','-r300');
savefig([filenameSIG '-2Dsurfplot.fig']);
clf();

%% make integrated plot 
%range for integration
range = length(x)/8;
subInt = Int(:,minX-range:minX+range);

intInt = sum(subInt,2);
if strcmp(subtract,'no')
    intInt = intInt - intInt(end); %subtract underground; this should be done with a underground measurement
end

if strcmp(plottype,'lin')
    plot(energy, intInt, 'linewidth',2);
    set(gca, 'Ylim',[0 max(intInt)+10000]);
    hold on;
%         %% fit the falling slope
%         [Max, maxIndex] = max(intInt);
%         timeFit = time(maxIndex:end)- time(maxIndex);
%         IntFit = intInt(maxIndex:end);
%         [f,gof,~] = fit(timeFit,IntFit,'exp1', 'StartPoint',[Max, -0.005]); 
%         %f(x) =  a*exp(b*x)
%         h = plot(timeFit+time(maxIndex), f.a*exp(f.b*timeFit));
%         h.LineWidth = 1.5;
%         h.Color = 'r' ; %[0 1 1]
%         l = legend('Data','Fit');
%         l.FontSize = 22;
%         decaytime = -1/f.b;
%         level = 2*tcdf(-1,gof.dfe);
%         m = confint(f,level); 
%         std = m(end,end) - f.b;
%         decaytimeError = 1/f.b^2 * std;
%         hold off;
%         set(gca,'DefaultTextInterpreter','latex');
%         text(0.4, 0.5, ['$\tau = $ ' num2str(decaytime,'%.1f') ' $\pm$ ' ...
%             num2str(decaytimeError,'%.1f') ' ps'],...
%             'FontSize',fontsize-4,...
%             'Units','normalized','Interpreter','latex');

elseif strcmp(plottype,'log')
    semilogy(energy(intInt>0), intInt(intInt>0), 'linewidth',2);
    hold on;
    %% fit the falling slope
%         [Max, maxIndex] = max(intInt);
%         IntFit0 = intInt(maxIndex:end);
%         IntFit = log(IntFit0(IntFit0>0));
%         range = max(IntFit)-min(IntFit);
%         level = range - 1;
%         levelIndex = find(max(IntFit)-IntFit>level,1,'first');
%         %IntFit = IntFit(1:levelIndex); %use only values before too much noise
%         timeFit = time(maxIndex:end)- time(maxIndex);
%         timeFit = timeFit(IntFit0>0);
%         %timeFit = timeFit(1:levelIndex);
%         [f,gof,~] = fit(timeFit,IntFit,'poly1'); %fity = p1*fitx + p2
%         b = f.p1; decaytime = -1/b;
%         level = 2*tcdf(-1,gof.dfe);
%         m = confint(f,level); 
%         std = m(end,1) - b;
%         decaytimeError = 1/b^2 * std;
%         h = plot(timeFit+time(maxIndex), exp(f.p2)*exp(b*timeFit));
%         h.LineWidth = 1.5;
%         h.Color = 'r' ; %[0 1 1]
%         l = legend('Data','Fit');
%         l.FontSize = 22;
%         hold off;
%         set(gca,'DefaultTextInterpreter','latex');
%         text(0.4,0.5, ['$\tau = $ ' num2str(decaytime,'%.1f') ' $\pm$ ' ...
%             num2str(decaytimeError,'%.1f') ' ps'],...
%             'FontSize',fontsize-4,...
%             'Units','normalized','Interpreter','latex');         
end
ylabel('Integrated Intensity (a.u.)','FontSize',fontsize,'FontName',fontName);
xlabel('Energy (eV)','FontSize',fontsize,'FontName',fontName);
set(gca,'LineWidth',3,'XColor',[0 0 0], 'YColor', [0 0 0],'Box','on',...
'FontSize',22,'FontName',fontName,...
'TickDir','In');
print([filenameSIG '-IntegratedPlot-' plottype '.png'],'-dpng','-r300');
savefig([filenameSIG '-IntegratedPlot-' plottype '.fig']);
clf();
    
end