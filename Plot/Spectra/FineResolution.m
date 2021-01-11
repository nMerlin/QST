%phasematching-fineresolution. Plots intensity values taken for different
%signal laser wavelenght written into an excel table. 

fitfunction = '';
%fitfunction = 'gauss';
%fitfunction = 'lorentz';

group = 16;

%overview
num = xlsread(['Phasematching-FineResolution-Group' num2str(group) '-2019-01-28.xlsx'],'A8:D130');
%num = xlsread(['Phasematching-FineResolution-Group' num2str(group) '-181212-results.xlsx'],'A8:D116');
wls =  num(:,1);
power =  num(:,2);
wlp = num(:,3);
%wli =  num(:,4);


%% plot 
% plot power
col = 'b';
hLine1 = plot(wls,power);
hLine1.LineStyle = 'none';
hLine1.Marker = 'o';
hLine1.MarkerEdgeColor = col;
hLine1.MarkerFaceColor = col;

hold on;

if strcmp(fitfunction,'gauss')
    [f,gof,~] = fit(wls,power,'gauss1', 'StartPoint',[60000, 897, 0.5]); %f(x) =  a1*exp(-((x-b1)/c1)^2)
    % f = fit(time,visibility,'gauss1');
    % gaussCustom = 'a1*exp(-((x-b1)/c1)^2)+d1';
    % f = fit( t, v, gaussCustom, 'StartPoint', [1, 0, 50,0.5] ); 
elseif strcmp(fitfunction,'lorentz')
    lorentz = 'a1 * c1/ (2*pi)/((x-b1)^2+(c1/2)^2)';
    [f,gof,~] = fit(wls,power,lorentz, 'StartPoint', [60000, 0, 0.5] );  
end

if strcmp(fitfunction,'gauss')||strcmp(fitfunction,'lorentz')
    FWHM = 2*f.c1*sqrt(log(2));
    cohTime = f.c1*sqrt(pi/2); 
    level = 2*tcdf(-1,gof.dfe);
    m = confint(f,level); 
    std = m(end,end) - f.c1;
    FWHMerror = std * 2*sqrt(log(2));
    cohTimeError = std*sqrt(pi/2); 
    
    h = plot(f);
    h.LineWidth = 1.5;
    h.Color = 'r' ; %[0 1 1]
    l = legend('Messdaten','Fit');
    l.FontSize = 22;
end

%interpolation
ip = csaps(wls, power);
fnplt(ip,'r-',0.5);

ylim([0 max(power)+500]);
fontsize = 24;
fontName = 'Times New Roman';
    set(gca,'LineWidth',3,'XColor',[0 0 0], 'YColor', [0 0 0],'Box','on',...
        'FontSize',22,'FontName',fontName,...
        'TickDir','Out');
ylabel('Converted Intensity (a.u.)','FontSize',fontsize,'Interpreter','latex')
xlabel('Signal (nm)','FontSize',fontsize,'Interpreter','latex') % label x-axis
legend('data','spline','Location','northwest');
savefig(['fineResolution-group' num2str(group) '-pump-' num2str(mean(wlp(~isnan(wlp)))) fitfunction '.fig']);
print(['fineResolution-group' num2str(group) '-pump-' num2str(mean(wlp(~isnan(wlp)))) fitfunction '.png'], '-dpng');
hold off;

