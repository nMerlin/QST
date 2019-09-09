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
[subtract] = c{:};

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
 
%% make plot
surf(k, energy, log(Int));
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
    
    
end