function [Max, integratedInt, peak, FWHM, Q, duration,FWHMdata] = plotSpectrumAndFit(filenameSIG,filenameBG,varargin)
%% PLOTSPECTRUM plots an optical spectrum, subtracts background, 
% fits a Gaussian to the peak. 
%
%   Input Arguments:
%       filenameSIG: file with the data. The file should be located in folder
%       'raw-data'. The data should consist of two columns; 
%       first column contains wavelength and second column contains
%       intensities
%       filenameBG: file with background data.
%       'Fit': Decide with 'yes', if the data should be
%       fitted with a gaussian.
%       'Interpolate': Decide with 'yes', if the data should be
%       interpolated with a spline. This interpolation will also be used
%       for the fit.
%       'XLim': limits for x-axis around the peak.
%       'XUnit': x unit can be nm (wavelength) or Hz (frequency).

%% Validate and parse input arguments
parser = inputParser;
defaultFit = 'yes'; %
addParameter(parser,'Fit',defaultFit);
defaultSave = 'yes'; %
addParameter(parser,'Save',defaultSave);
defaultInterpolate = 'yes';
addParameter(parser,'Interpolate',defaultInterpolate);
defaultSubtract = 'no'; %
addParameter(parser,'Subtract',defaultSubtract);
defaultXLim = 0.5; %
addParameter(parser,'XLim',defaultXLim,@isnumeric);
defaultXUnit = 'nm'; % alternatives: Hz, eV
addParameter(parser,'XUnit',defaultXUnit);
defaultXrange = [];
addParameter(parser,'Xrange',defaultXrange);
defaultNormalize = false;
addParameter(parser,'Normalize',defaultNormalize,@islogical);
defaultAnnotate = true;
addParameter(parser,'Annotate',defaultAnnotate,@islogical);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[annotate,fitoption,intp,normalize,save,subtract,xLim,xRange,xUnit] = c{:};

%%
fontsize = 18;
lightVelocity = 299792458;
h = 6.62607015e-34;
e0 = 1.602176634e-19;

%% load data
    cd('raw-data')
   data = textread(filenameSIG);
   % data = textread(filenameSIG,'','delimiter',','); 
  % data = textread(filenameSIG,'','delimiter',',','headerlines',1);
  
%for data from stellarnet spectrometer:
%     x = importdata(filenameSIG,' ',2);
%     data = x.data;  
    
    w = data(:,1); % wavelength
    Int = data(:,2); % Intensity
    %Int = data(:,3); % Intensity --> sometimes the intensity is in the 3rd column 
    
%% confine data in wavelength range if wished
if not(isempty(xRange))
    Int = Int((w>= min(xRange)) & (w<= max(xRange))); 
    w = w((w>= min(xRange)) & (w<= max(xRange)));
end
    
%% subtract background
    if strcmp(subtract, 'yes')
        dataBG = textread(filenameBG);
        IntBG = dataBG(:,2);
        Int = Int - IntBG;
    else
%         bg = min(Int);
%         Int = Int-bg;
    end

%% normalize data
if normalize
   Int = Int ./max(Int); 
end
    
%% get maximum and peak position
    [Max,I] = max(Int);
    integratedInt = sum(Int);
    peak = w(I);
    
    if strcmp(xUnit,'Hz')
        xLim = lightVelocity*xLim*10^-9./(peak*10^-9)^2;
        w = lightVelocity./(w*10^-9);
        peak = w(I);
    elseif strcmp(xUnit,'eV')
        xLim = h*lightVelocity/e0*xLim*10^-9./(peak*10^-9)^2;
        w = h*lightVelocity/e0*1./(w*10^-9);
        peak = w(I);
    end
    

    
%% interpolate data, if wished
    if strcmp(intp, 'yes')
        xx= min(w):xLim/1000:max(w);%xx= peak-xLim:0.001:peak+xLim;
        ip = spline(w, Int, xx);
        Intfit = ip(xx>peak-xLim & xx<peak+xLim)';
        wfit = xx(xx>peak-xLim & xx<peak+xLim)';
        plot(xx,ip,'LineWidth',1,'DisplayName','Interp.'); hold on;
        plot(w,Int,'bo','LineWidth',2,'DisplayName','Data');
        hold on;
    else
        Intfit = Int;
        wfit = w;
        plot(w,Int,'b-','LineWidth',2,'DisplayName','Data');
        hold on;
    end    
    
%% FWHM from data  
    FWHMData = fwhm(w,Int);
        if strcmp(xUnit, 'nm')
            durationData= 2*log(2)/pi *(peak*1e-9)^2 / lightVelocity / (FWHMData*1e-9) * 1e15; %time in femtoseconds
        elseif strcmp(xUnit, 'eV')
            durationData = 2*log(2)/pi * h/e0/FWHMData *1e15;
        elseif strcmp(xUnit, 'Hz')
            durationData = 2*log(2)/pi / FWHMData *1e15;
        end 
        if annotate
            text(peak-0.9*xLim, Max,...
            ['FWHMData = ' num2str(FWHMData,'%.4f') ' ' xUnit char(10) ...
            '\Delta t = ' num2str(durationData,'%.1f') ' fs'],'FontSize',fontsize-4);
        end
    
%% make Gauss Fit if wished
    if strcmp(fitoption, 'yes')
        %gaussCustom = 'a1*exp(-((x-b1)/c1)^2)+d1';
        [f,gof,~] = fit(wfit,Intfit,'gauss1', 'StartPoint', [Max, peak, xLim] ); %f(x) =  a1*exp(-((x-b1)/c1)^2)
        FWHM = 2*f.c1*sqrt(log(2));
        % Properties of gaussian: https://de.wikipedia.org/wiki/Normalverteilung
        level = 2*tcdf(-1,gof.dfe);
        m = confint(f,1-level); 
        std = m(end,end) - f.c1;
        FWHMerror = std * 2*sqrt(log(2));
        %relative width
        Q = f.b1/FWHM;
        % In order to compute the duration in time, the duration in
        % whatever unit is converted to a duration in frequency space and
        % then we use the relation between duration in frequency and time
        % domain given by Fourier transformation.
        % According to
        % https://de.wikipedia.org/wiki/Fourier-Transformation#Tabelle_wichtiger_Fourier-Transformations-Paare,
        % we get sigma_t = 1/(2*pi*sigma_f). Also we use the relation
        % sigma = FWHM / (2*sqrt(log(2)). 
        % Then : FWHM_t = 2*log(2)/pi / FWHM_f.
        % Another source for this: https://de.wikipedia.org/wiki/Pulslaser
        if strcmp(xUnit, 'nm')
            duration= 2*log(2)/pi *(peak*1e-9)^2 / lightVelocity / (FWHM*1e-9) * 1e15; %time in femtoseconds
        elseif strcmp(xUnit, 'eV')
            duration = 2*log(2)/pi * h/e0/FWHM *1e15;
        elseif strcmp(xUnit, 'Hz')
            duration = 2*log(2)/pi / FWHM *1e15;
        end               
        x= peak-xLim:xLim/1000:peak+xLim; %x = peak-xLim:0.001:peak+xLim;
        plot(x,f(x),'r','LineWidth',1.5,'DisplayName','Fit');
        %set(gca,'DefaultTextInterpreter','latex');
        if annotate
            text(peak-0.9*xLim, Max/4,...
            ['FWHMGauss = ' num2str(FWHM,'%.4f') ' \pm ' num2str(FWHMerror,'%.4f') ' ' xUnit char(10) ...
            'Q = ' num2str(Q,'%.0f') char(10) '\Delta t = ' num2str(duration,'%.1f') ' fs'],'FontSize',fontsize-4);
        end
    else
        FWHM = 0;
        Q = 0;
    end
    
%% plot
    l = legend('Location','northeast');
    l.FontSize = 22;
    ylim([0 Max*1.1]);
    if not(isempty(xRange))
       xlim([min(xRange) max(xRange)]);
    else
        xlim([peak-xLim, peak+xLim]);
    end
    if normalize
        ylabel('Normalized Counts');
    else
        ylabel('Intensity (a.u.)');
    end    
    if strcmp(xUnit, 'nm')
        xlabel('Wavelength (nm)');    
    elseif strcmp(xUnit, 'Hz')
        xlabel('Frequency (Hz)');  
    elseif strcmp(xUnit, 'eV')
        xlabel('Energy (eV)');  
    end
    
    %set(gca,'DefaultTextInterpreter','latex');
    if annotate
        text(peak-0.9*xLim, Max/2,['peak at ' num2str(peak,'%.4f') ' ' xUnit char(10) ...
        'max Int ' num2str(Max,'%.0f') ' counts' char(10) ...
        'integr. Int ' num2str(integratedInt,'%.0f') ' counts' char(10)],...
        'FontSize',fontsize-4);
    end
    graphicsSettings('Fontsize',25);
    grid;
    hold off;
    cd('..')
    if strcmp(save, 'yes')
         savefig([filenameSIG '-plot-Subtract-' subtract '-Normalize-' num2str(normalize) '-annotate-' num2str(annotate) '-' xUnit '.fig']);
         print([filenameSIG '-plot-Subtract-' subtract '-Normalize-' num2str(normalize) '-annotate-' num2str(annotate) '-' xUnit '.png'],'-dpng','-r300');
    end

end
