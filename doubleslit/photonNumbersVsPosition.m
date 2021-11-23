function [fittedPeriodInMM,dInM] = photonNumbersVsPosition(varargin)

%% Validate and parse input arguments
p = inputParser;
defaultFolder = 'mat-data'; % which folder contains the data
addParameter(p,'Folder',defaultFolder);
parse(p,varargin{:});
c = struct2cell(p.Results);
[folder] = c{:};


[filenames,numbers,Positions]= getParametersFromFilenames('Folder',folder,'Parameter','position');
photonnumbers = zeros(length(filenames),1);
for fileI = 1:length(filenames)
%     load([folder '\' cell2mat(filenames(fileI))],'X2');
%     [~,n,~] = nPhotons(X2,X2,X2);
    load([folder '\' cell2mat(filenames(fileI))],'X1');
    [~,n,~] = nPhotons(X1,X1,X1);
    X2 = X1;
    photonnumbers(fileI) = n;
    X2 = X2(:);
%     plot(X2(:),'.');
%     savefig([cell2mat(filenames(fileI)) '-X.fig']);
%     clf();
end

numberOfPeriods = 9;
plot(Positions,photonnumbers,'o-');
xlabel('Position (mm)');
ylabel('Mean photon number');
graphicsSettings();
savefig('nVsPosition.fig');
print('nVsPosition.png','-dpng','-r300');
clf();

 %% normalize the interference pattern to absolute maximum and minimum 
y = photonnumbers;
[ynew,locs,maxlocs,minlocs] = normalizeInterferencePattern(y,numberOfPeriods);
plot(Positions,ynew,'o');
graphicsSettings();
hold on;
plot(Positions,ynew,'b-');
xlabel('Position (mm)');
ylabel('normalized mean photon number');
%plot theory curve 
%% Estimate fit paramters
y = ynew;
yMax = mean(y(maxlocs));
yMin = mean(y(minlocs));
yRange = (yMax-yMin);
xphase = -Positions(minlocs(1));
period  = mean(diff(Positions(locs)))*2;
fitFunction = 'a.*(sin(2*pi*x./b + c)).^2';
[f,gof,~] = fit(Positions',y,fitFunction, 'StartPoint', [yRange, period*2, xphase] );
hold on;
xNew = 0:0.01:Positions(end);
plot(xNew,yRange.*(sin(2*pi*xNew./(f.b) + f.c)).^2 + yMin,'r-','LineWidth',2);
fittedPeriodInMM = f.b/2;
savefig('nVsPosition-normalized.fig');
print('nVsPosition-normalized.png','-dpng','-r300');
clf();

%Interferenzwegl?nge in m
D = 0.6;
%Wellenl?nge in m
lambda = 834e-9;
%Abstand Maxima in m 
a = fittedPeriodInMM*1e-3;
%Abstand der Spalte
dInM = D*lambda/a;

save('photonNumbersVsPosition.mat','filenames','numbers','delays','photonnumbers','fittedPeriodInMM','dInM');

end 