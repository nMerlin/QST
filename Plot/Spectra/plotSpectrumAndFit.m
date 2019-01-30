function [Max, peak, FWHM, Q] = plotSpectrumAndFit(filename,varargin)
%% PLOTSPECTRUM plots an optical spectrum, subtracts background, 
% fits a Gaussian to the peak. 
%
%   Input Arguments:
%       filename: file with the data. The file should be located in folder
%       'raw-data'. The data should consist of two columns; 
%       first column contains wavelength and second column contains
%       intensities
%       'Interpolate': Decide with 'yes', if the data should be
%       interpolated with a spline. This interpolation will also be used
%       for the fit.

%% Validate and parse input arguments
parser = inputParser;
defaultInterpolate = 'yes'; %
addParameter(parser,'Interpolate',defaultInterpolate);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[intp] = c{:};

%% load data
    cd('raw-data')
    data = textread(filename);
    
    w = data(:,1); % wavelength
    Int = data(:,2); % Intensity
    
%% subtract background
    bg = min(Int);
    Int = Int-bg;
    
%% get maximum and peak position
    [Max,I] = max(Int);
    integratedInt = sum(Int);
    peak = w(I);
    
%% interpolate data, if wished
    if strcmp(intp, 'yes')
        xx= peak-2:0.001:peak+2;
        ip = spline(w, Int, xx);
        Intfit = ip';
        wfit = xx';
        plot(xx,ip,'LineWidth',1); hold on;
    else
        Intfit = Int;
        wfit = w;
    end    
    
%% make Gauss Fit
    %gaussCustom = 'a1*exp(-((x-b1)/c1)^2)+d1';
    [f,gof,~] = fit(wfit,Intfit,'gauss1', 'StartPoint', [Max, peak, 0.5] ); %f(x) =  a1*exp(-((x-b1)/c1)^2)
    FWHM = 2*f.c1*sqrt(log(2));
    level = 2*tcdf(-1,gof.dfe);
    m = confint(f,level); 
    std = m(end,end) - f.c1;
    FWHMerror = std * 2*sqrt(log(2));
    %relative width
    Q = f.b1/FWHM;
    
%% plot
    plot(w,Int,'bo','LineWidth',2);
    hold on;
    x = peak-2:0.001:peak+2;
    plot(x,f(x),'r','LineWidth',1.5);
    l = legend('Interp.','Data','Fit');
    l.FontSize = 22;
    ylim([0 Max*1.1]);
    xlim([peak-0.5, peak+0.5]);
    ylabel('Counts');
    xlabel('Wavelength (nm)');
    fontsize = 24;
    fontName = 'Times New Roman';
    set(gca,'LineWidth',3,'XColor',[0 0 0], 'YColor', [0 0 0],'Box','on',...
        'FontSize',22,'FontName',fontName,...
        'TickDir','Out');
    set(gca,'DefaultTextInterpreter','latex');
    text(peak-0.45, Max/2,['peak at ' num2str(peak,'%.2f') ' nm' char(10) ...
        'max Int ' num2str(Max,'%.0f') ' counts' char(10) ...
        'integr. Int ' num2str(integratedInt,'%.0f') ' counts' char(10) ...
        'FWHM = ' num2str(FWHM,'%.3f') ' $\pm$ ' num2str(FWHMerror,'%.4f') ' nm' char(10) ...
        'Q = ' num2str(Q,'%.0f')],...
        'FontSize',fontsize-4);
    hold off;
    cd('..')
    savefig([filename '-subtractedBackground-plot.fig']);
    print([filename '-subtractedBackground-plot.png'],'-dpng','-r300');

end
