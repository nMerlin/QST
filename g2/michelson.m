function [ t,v, FWHM, FWHMerror, cohTime, cohTimeError] = michelson(current,varargin)
%MICHELSON Creates a plot showing the visibility over time delay of a
%michelson interferometer measurement.
%
% Input Arguments:
%   current: current in mA. Helps identify the right overview file.
%   This function needs an overview file called 'currentmA-overview.xlsx' 
%   assigning a delay line position in counts to each oscilloscope measurement.
%   The oscilloscope filed have to be in a folder 'raw-data'.
%
% Optional Input Arguments:
%   'Fit': the fitfunction. Can be 'gauss', 'lorentz' or none. 
%
% Output Arguments:
%   t: time delay
%   v: visibility of the interference
%   Unit: units for the time can be nano 'n' or femto 'f' seconds.
%   FWHM, FWHMerror, cohTime, cohTimeError: FWHM and its fitting error. Coherence
%   time and its fitting error.

%% Validate and parse input arguments
p = inputParser;
defaultFit = 'gauss';
addParameter(p,'Fit',defaultFit);
defaultUnit = 'f';
addParameter(p,'Unit',defaultUnit);
parse(p,varargin{:});
c = struct2cell(p.Results);
[fitfunction, unit] = c{:};

%overview
[position,filenames] = xlsread([num2str(current) 'mA-overview.xlsx'],'A2:B100');

%position from counts to m
positionM = position / 2048000;
%time from m to seconds
c = 299792458; % in m/s
t = 2* positionM / c; %factor 2 because the beam goes the way two times
 %converts into wished unit
if strcmp(unit,'n')
    t = t * 10^9;
elseif strcmp(unit,'f')
    t = t * 10^15;
end

%% compute visibiltiy from oscilloscope traces 
v = zeros(length(position),1);
cd('raw-data');

for i = 1:length(position)

    m = csvread(cell2mat(filenames(i)),18); %old oscilloscope
    P = m(:,5); %old oscilloscope %measured signal
%     m = csvread(cell2mat(filenames(i)),2); %new oscilloscope
%     P = m(:,1); %new oscilloscope
    x = 1:length(P);
    Ps = csaps(x,P,1e-4,x); %smoothing
    [maxpks, ~] = findpeaks(Ps,'minPeakProminence',0.02*abs(max(Ps)));
    [minpks, ~] = findpeaks(-Ps,'minPeakProminence',0.02*abs(max(-Ps)));
    minpks = -minpks;
    findpeaks(Ps,'minPeakProminence',0.01*abs(max(Ps)));
    findpeaks(-Ps,'minPeakProminence',0.01*abs(max(-Ps)));
    
    %take highest maximum and smallest minimum
    maxP = max(maxpks);
    minP = min(minpks);
    
    if isempty(maxP) || isempty(minP)
        v(i) = 0;
    else
        v(i) = (maxP-minP)/(maxP+minP);
    end
end

cd('..');
%t = t-mean(t); % when the position with max. visibility is not zero.

%sort
[t,I] = sort(t);
v = v(I);

%% Plotting and fitting
plot(t,v,'ok','linewidth',1.5); %'ob'
hold on;

if strcmp(fitfunction,'gauss')
    [f,gof,~] = fit(t,v,'gauss1', 'StartPoint', [1, 0, 100] ); %f(x) =  a1*exp(-((x-b1)/c1)^2)
    % f = fit(time,visibility,'gauss1');
    % gaussCustom = 'a1*exp(-((x-b1)/c1)^2)+d1';
    % f = fit( t, v, gaussCustom, 'StartPoint', [1, 0, 50,0.5] ); 
elseif strcmp(fit,'lorentz')
    lorentz = 'a1 * c1/ (2*pi)/((x-b1)^2+(c1/2)^2)';
    [f,gof,~] = fit(t,v,lorentz, 'StartPoint', [1, 0, 100] );  
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

ylim([0 1]);
fontname = 'Times New Roman';

%title(['I = ' num2str(current) ' mA; {\tau_{c}} = ' num2str(cohTime) ' fs']);
%title(['I = ' num2str(current) ' mA; FWHM = ' num2str(FWHM) ' fs']);
% title(['I = ' num2str(current) ' mA']);
fontsize = 24;
xlabel(['$\tau$ (' unit 's)'],'FontSize',fontsize,'Interpreter','latex');
ylabel('$ g^{(1)}( \tau) $','FontSize',fontsize,'Interpreter','latex');

set(gca, 'LineWidth',1.5, 'FontSize', 22, 'FontName',fontname, 'XColor', 'k', 'YColor', 'k'); 

print(['Visibility-' num2str(current) 'mA-MA.png'], '-dpng');
% bw = imread(strcat('Visibility-',num2str(current),'mA-FWHM-',num2str(FWHM),'fs.png'));
% bw2 = imcomplement(bw);
% imwrite(bw2,strcat('Visibility-',num2str(current),'mA-FWHM-',num2str(FWHM),'fs-inverted.png'));
    
end
