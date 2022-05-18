function [PthFromModeInt,PthErrFromModeInt,PthFromFit,PthErrFromFit] = plotDispersionsPowerSeries(varargin)
% this script makes spectrum plots for a series of dispersion measurements.
% The data should be in folder 'raw-data'.
% With 'Subtract', you can choose whether there are background measurements
% that should be subtracted. 

%% Validate and parse input arguments
parser = inputParser;
defaultSubtract = 'no'; %
addParameter(parser,'Subtract',defaultSubtract);
defaultPlottype = 'lin'; 
addParameter(parser,'Plottype',defaultPlottype);
defaultPlotMode = 'no'; 
addParameter(parser,'PlotMode',defaultPlotMode); % if yes, the selected mode is indicated in the plot. 
defaultMinY = 200; 
addParameter(parser,'minY',defaultMinY,@isnumeric);
defaultXAperture = 330'; 
addParameter(parser,'xAperture',defaultXAperture,@isnumeric);
defaultModeK = [-0.66 0.66]; 
addParameter(parser,'ModeK',defaultModeK,@isnumeric); % k range in which the mode is selected
defaultModeE = [1.611 1.6113]; 
addParameter(parser,'ModeE',defaultModeE,@isnumeric); % Energy range in which the mode is selected
defaultAdjustEnergy = false; %set whether the energy range of the mode is adjusted to lie around the fitted peak for computing the integrated mode intensity 
addParameter(parser,'AdjustEnergy',defaultAdjustEnergy);
defaultZoomE = [1.605 1.625]; 
addParameter(parser,'ZoomE',defaultZoomE,@isnumeric); % Energy range in which the dispersion is plotted
defaultOnlySmallRange = false;
addParameter(parser,'OnlySmallRange',defaultOnlySmallRange,@islogical); %get FWHM and peakPosition only in the selected energy range
defaultFit = 'yes'; %'useOld' to use old parameters
addParameter(parser,'Fit',defaultFit);
defaultE0Old = 1.6108; 
addParameter(parser,'E0Old',defaultE0Old,@isnumeric);
defaultaOld = 4.4973e-7; 
addParameter(parser,'aOld',defaultaOld,@isnumeric);
defaultY0Old = 217; 
addParameter(parser,'y0Old',defaultY0Old,@isnumeric);
defaultGetTime =false; %
addParameter(parser,'GetTime',defaultGetTime,@islogical); % get exposure time from filename
defaultLabel = '-'; %
addParameter(parser,'Label',defaultLabel); % use only files with specific label e.g. polarisation
defaultGetThresholdFromFit = false; %get the threshold power from fit
addParameter(parser,'GetThresholdFromFit',defaultGetThresholdFromFit,@islogical);
defaultPlotOption = true;
addParameter(parser,'PlotOption',defaultPlotOption,@islogical);
defaultNormalize = false; %normalizes the colorscale in the 2D plot to 1.
addParameter(parser,'Normalize',defaultNormalize,@islogical);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[adjustEnergy,aOld,E0Old,fitoption,getThresholdFromFit,getTime,label,minY,modeE,modeK,...
    normalize,onlySmallRange,plotMode,plotOption,plottype,subtract,xAperture,y0Old,zoomE] = c{:};
%for filename:
options = ['-getTime-' num2str(getTime) '-modeE-' num2str(min(modeE)) '-' num2str(max(modeE)) ...
     '-modeK-' num2str(min(modeK)) '-' num2str(max(modeK)) '-adjustModeEnergy-' num2str(adjustEnergy),'-onlySmallRange-',num2str(onlySmallRange)];

%% Sample char. for M3396
R = 9.5; %Rabi splitting in meV
EX = 1.6195; % Exciton energy in eV
%% constants
h = 6.62607015e-34;
c0 = 299792458;
e0 = 1.602176634e-19;

%% Create data overview
dataStruct = struct('filename',{},'number',{}, 'Power',{}, 'time', {}, 'E0', {}, 'lambda',{}, 'Emax', {}, 'modeInt', {}, 'SumInt', {});
dataStructBackground = struct('filename',{},'number',{});

rawDataContents = dir('raw-data');

%% Signal-files
for name = {rawDataContents.name}
    % Loop only over Signal files
    filename = cell2mat(name);
    if not(isempty(regexpi(filename,'-background.csv','match')))...
            || not(isempty(regexpi(filename,'raw','match')))...
            || isempty(regexpi(filename,'.csv','match'))...
            || isempty(regexpi(filename,label,'match'))
        continue
    end
    
    % Fetch number of measurement
    numberToken = regexpi(filename,'^([0123456789]*)-','tokens');
    number = str2double(cell2mat(numberToken{1}));
    dataStruct(number).filename = filename;
    
    % Fetch excitation power
    powerToken = regexpi(filename,'-([0123456789.]*)mW','tokens');
    power = str2double(cell2mat(powerToken{1}));
    dataStruct(number).Power = power;
    
    if getTime
        timeToken = regexpi(filename,'-([0123456789.]*)ms','tokens');
        time = str2double(cell2mat(timeToken{1}));
        dataStruct(number).time = time;
    else
        dataStruct(number).time = 1;
    end
    
end

%% Background-files
for name = {rawDataContents.name}
    % Loop only over background files
    filename = cell2mat(name);
    if isempty(regexpi(filename,'background','match'))
        continue
    end
    
    % Fetch number of measurement
    numberToken = regexpi(filename,'^([0123456789]*)-','tokens');
    number = str2double(cell2mat(numberToken{1}));
    dataStructBackground(number).filename = filename;
    dataStructBackground(number).number = number;
end
BGnumbers = cell2mat({dataStructBackground.number});

%% process the data
for number = 1:size(dataStruct,2)

    filenameSIG = dataStruct(number).filename;
    if isempty(filenameSIG)
        continue
    end    
   
    %find adequate BG-file, if this option is chosen
    if strcmp(subtract, 'yes')
        BGnumber = min(BGnumbers(BGnumbers<=number)); %background was measured before signal
        filenameBG = dataStructBackground(BGnumber).filename;
    else
        filenameBG = filenameSIG;
    end
    [E0,~,~,Emax,modeInt,SumInt,FWHMFit,FWHMerror,peakPositionFit,peakHeight,~,peakPosition,FWHM,Areapixels] = plotDispersion(filenameSIG, filenameBG,'Subtract',subtract,...
        'Plottype',plottype,'minY',minY,'xAperture',xAperture,'ModeE',modeE,...
        'ModeK',modeK,'ZoomE',zoomE,'Fit',fitoption,'aOld',aOld,'E0Old',E0Old,...
        'y0Old',y0Old,'PlotMode',plotMode,'R',R,'EX',EX,'PlotOption',plotOption,...
        'Normalize',normalize,'AdjustEnergy',adjustEnergy,'OnlySmallRange',onlySmallRange);
    lambda = h*c0/e0*1./E0 *10^9;
    dataStruct(number).E0 = E0;
    dataStruct(number).lambda = lambda;
    dataStruct(number).Emax = Emax; 
    dataStruct(number).FWHMFit =FWHMFit; 
    dataStruct(number).FWHMerror = FWHMerror; 
    dataStruct(number).FWHM =FWHM; 
    dataStruct(number).peakPositionFit = peakPositionFit; 
    dataStruct(number).peakPosition = peakPosition; 
    dataStruct(number).Areapixels = Areapixels; 
    if getTime
        dataStruct(number).modeInt = modeInt/dataStruct(number).time;
        dataStruct(number).SumInt = SumInt/dataStruct(number).time;
        dataStruct(number).peakHeight = peakHeight/dataStruct(number).time;
    else
        dataStruct(number).modeInt = modeInt; 
        dataStruct(number).SumInt = SumInt;
        dataStruct(number).peakHeight = peakHeight;
    end
    
end

power = cell2mat({dataStruct.Power});
E0 = cell2mat({dataStruct.E0});
lambda = cell2mat({dataStruct.lambda});
Emax = cell2mat({dataStruct.Emax});
modeInt = cell2mat({dataStruct.modeInt});
SumInt = cell2mat({dataStruct.SumInt});
FWHMFit = cell2mat({dataStruct.FWHMFit});
FWHMerror = cell2mat({dataStruct.FWHMerror});
FWHM = cell2mat({dataStruct.FWHM});
peakPositionFit = cell2mat({dataStruct.peakPositionFit});
peakPosition = cell2mat({dataStruct.peakPosition});
peakHeight = cell2mat({dataStruct.peakHeight});
Areapixels = cell2mat({dataStruct.Areapixels});
time = cell2mat({dataStruct.time});
save(['powerSeries' label options '.mat'],'power','E0','lambda','Emax','modeInt','SumInt','FWHMFit',...
    'FWHMerror','peakPositionFit','peakHeight','time','peakPosition','FWHM','Areapixels');

%% write them in excel table
T = struct2table(dataStruct);
writetable(T,['powerSeries' label options '.xls']);

%% make plots
fontsize = 20;
fontname = 'Arial';

%%
semilogx(power,peakPositionFit,'ko','markerFaceColor','k');
xlabel('Excitation Power P (mW)');
ylabel('Energy (eV)');
graphicsSettings;
ax = gca;
set(ax,'FontSize',fontsize,'FontName',fontname);
savefig(['powerSeries-peakPositionFit' label options '.fig']);
print(['powerSeries-peakPositionFit' label options '.png'],'-dpng','-r200');
clf();

%%
semilogx(power,peakPosition,'ko','markerFaceColor','k');
xlabel('Excitation Power P (mW)');
ylabel('Energy (eV)');
graphicsSettings;
ax = gca;
set(ax,'FontSize',fontsize,'FontName',fontname);
savefig(['powerSeries-peakPosition' label options '.fig']);
print(['powerSeries-peakPosition' label options '.png'],'-dpng','-r200');
clf();

%%
loglog(power,peakHeight,'ko','markerFaceColor','k');
xlabel('Excitation Power P (mW)');
if getTime
   ylabel('Peak Height of Fit (counts/ms)'); 
else
    ylabel('Peak Height of Fit (counts)');
end
graphicsSettings;
ax = gca;
set(ax,'FontSize',fontsize,'FontName',fontname);
if getThresholdFromFit
    xfit1 = log(power(1:5));
    yfit1 = log(peakHeight(1:5));
    fitFunction = fittype('a*x + b');  
    [f,gof,~] = fit(xfit1',yfit1',fitFunction);
    hold on;
    loglog(power,exp(f(log(power))),'r-','Linewidth',2);
    a1 = f.a;    
    b1 = f.b;
    [se]= getStandardErrorsFromFit(f,gof,'method1');   
    a1Err = se(1);
    b1Err = se(2);
    xfit2 = log(power(6:9));
    yfit2 = log(peakHeight(6:9)); 
    [f2,gof2,~] = fit(xfit2',yfit2',fitFunction);
    hold on;
    loglog(power,exp(f2(log(power))),'r--','Linewidth',2);
    ylim([1e-2 1e6]);
    a2 = f2.a;    
    b2 = f2.b;
    [se2]= getStandardErrorsFromFit(f2,gof2,'method1');   
    a2Err = se2(1);
    b2Err = se2(2);
    PthFromFit = exp((b2-b1)/(a1-a2));
    %get error with montecarlo
    PthRand = zeros(1000,1);
    for i = 1:1000
        a1r = normrnd(a1,a1Err);
        a2r = normrnd(a2,a2Err);
        b1r = normrnd(b1,b1Err);
        b2r = normrnd(b2,b2Err);
        PthRand(i) = exp((b2r-b1r)/(a1r-a2r));
    end
    PthErrFromFit = std(PthRand);
     xfit3 = log(power(10:end));
    yfit3 = log(peakHeight(10:end)); 
    [f3,gof3,~] = fit(xfit3',yfit3',fitFunction);
    hold on;
    loglog(power,exp(f3(log(power))),'r-.','Linewidth',2);
    ylim([1e-2 1e6]);
    a3 = f3.a;    
    b3 = f3.b;
    [se3]= getStandardErrorsFromFit(f3,gof3,'method1');   
    a3Err = se3(1);
    b3Err = se3(2);
else
    PthFromFit = 0;
    PthErrFromFit = 0;
end
savefig(['powerSeries-peakHeight' label options '.fig']);
print(['powerSeries-peakHeight' label options '.png'],'-dpng','-r200');
clf();

%%
loglog(power,modeInt,'ko','markerFaceColor','k');
xlabel('Excitation Power P (mW)');
if getTime
    ylabel('Integrated Mode Intensity (counts/ms)');
else
    ylabel('Integrated Mode Intensity (counts)');
end
graphicsSettings;
ax = gca;
set(ax,'FontSize',fontsize,'FontName',fontname);
if getThresholdFromFit
    xfit1 = log(power(1:5)); %adjust these!
    yfit1 = log(modeInt(1:5));
    fitFunction = fittype('a*x + b');  
    [f,gof,~] = fit(xfit1',yfit1',fitFunction);
    hold on;
    loglog(power,exp(f(log(power))),'r-','Linewidth',2);
    a1 = f.a;    
    b1 = f.b;
    [se]= getStandardErrorsFromFit(f,gof,'method1');   
    a1Err = se(1);
    b1Err = se(2);
    xfit2 = log(power(6:9));
    yfit2 = log(modeInt(6:9)); 
    [f2,gof2,~] = fit(xfit2',yfit2',fitFunction);
    hold on;
    loglog(power,exp(f2(log(power))),'r--','Linewidth',2);
    a2 = f2.a;    
    b2 = f2.b;
    [se2]= getStandardErrorsFromFit(f2,gof2,'method1');   
    a2Err = se2(1);
    b2Err = se2(2);
    PthFromModeInt = exp((b2-b1)/(a1-a2));
    %get error with montecarlo
    PthRand = zeros(1000,1);
    for i = 1:1000
        a1r = normrnd(a1,a1Err);
        a2r = normrnd(a2,a2Err);
        b1r = normrnd(b1,b1Err);
        b2r = normrnd(b2,b2Err);
        PthRand(i) = exp((b2r-b1r)/(a1r-a2r));
    end
    PthErrFromModeInt = std(PthRand);
        ylim([1e-2 1e6]);
else
    PthFromModeInt = 0;
    PthErrFromModeInt = 0;
end
savefig(['powerSeries-modeInt' label options '.fig']);
print(['powerSeries-modeInt' label options '.png'],'-dpng','-r300');
clf();

%% modeInt normalized to number of pixels 
loglog(power,modeInt./Areapixels,'o-');
xlabel('Excitation Power P (mW)');
if getTime
    ylabel('Mode Intensity (counts ms^{-1} px^{-1})');
else
    ylabel('Mode Intensity (counts/px)');
end
graphicsSettings;
ax = gca;
set(ax,'FontSize',fontsize,'FontName',fontname);
savefig(['powerSeries-modeIntPerPixel' label options '.fig']);
print(['powerSeries-modeIntPerPixel' label options '.png'],'-dpng','-r300');
clf();

%%
semilogx(power,E0,'o');
xlabel('Excitation Power (mW)');
ylabel('minimum energy E_{LP}(k_{||}=0)(eV)');
graphicsSettings;
savefig(['powerSeries-E0' label options '.fig']);
print(['powerSeries-E0' label options '.png'],'-dpng','-r300');
clf();

semilogx(power,Emax,'o');
xlabel('Excitation Power (mW)');
ylabel('Energy of maximum of mode (eV)');
graphicsSettings;
savefig(['powerSeries-Emax' label options '.fig']);
print(['powerSeries-Emax' label options '.png'],'-dpng','-r300');
clf();

semilogx(power,lambda,'o');
xlabel('Excitation Power (mW)');
ylabel('wavelength of maximum of mode (nm)');
graphicsSettings;
savefig(['powerSeries-wavelength' label options '.fig']);
print(['powerSeries-wavelength' label options '.png'],'-dpng','-r300');
clf();

errorbar(power,FWHMFit,FWHMerror);
f = gca;
f.XScale = 'log';
xlabel('Excitation Power (mW)');
ylabel('FWHM of Fit around k = 0 (meV)');
graphicsSettings;
savefig(['powerSeries-FWHMFit' label options '.fig']);
print(['powerSeries-FWHMFit' label options '.png'],'-dpng','-r300');
clf();

semilogx(power,FWHM,'-o');
xlabel('Excitation Power (mW)');
ylabel('FWHM (meV)');
graphicsSettings;
savefig(['powerSeries-FWHM' label options '.fig']);
print(['powerSeries-FWHM' label options '.png'],'-dpng','-r300');
clf();

loglog(power,SumInt,'o');
xlabel('Excitation Power (mW)');
if getTime
   ylabel('overall integrated Intensity (counts/ms)'); 
else
    ylabel('overall integrated Intensity (counts)');
end
graphicsSettings;
savefig(['powerSeries-SumInt' label options '.fig']);
print(['powerSeries-SumInt' label options '.png'],'-dpng','-r300');
clf();


end