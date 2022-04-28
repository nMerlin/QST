function [] = correctedPlot(label)

% get data
% M = xlsread('powerSeries-filterCorrected.xls','A2:L50');
M = xlsread(['powerSeries' label 'filterCorrected.xls'],'A2:L50');
% power = M(:,6);
% Intmax = M(:,7);
% SumInt = M(:,8);
power = M(:,8);
modeInt = M(:,9);
SumInt = M(:,10);
%remove NaNvalues
power = power(~isnan(power));
modeInt = modeInt(~isnan(modeInt));
SumInt = SumInt(~isnan(SumInt));

loglog(power,modeInt,'o','Linewidth',2);
xlabel('Excitation Power (mW)');
%ylabel('maximum Intensity at k = 0 (a.u.)');
ylabel('mode Intensity at k = 0 (counts/time)');
graphicsSettings;
savefig(['powerSeries-modeInt' label 'filterCorrected.fig']);
print(['powerSeries-modeInt' label 'filterCorrected.png'],'-dpng','-r300');
clf();

% plot(power,Intmax,'o');
% xlabel('Excitation Power (mW)');
% %ylabel('maximum Intensity at k = 0 (a.u.)');
% ylabel('maximum Intensity at k = 0 (a.u.)');
% graphicsSettings;
% savefig('powerSeries-Intmax-filterCorrected-lin.fig');
% print('powerSeries-Intmax-filterCorrected-lin.png','-dpng','-r300');
% clf();

loglog(power,SumInt,'o','Linewidth',2);
xlabel('Excitation Power (mW)');
%ylabel('overall integrated Intensity (a.u.)');
ylabel('overall integrated Intensity (counts/time)');
graphicsSettings;
savefig(['powerSeries-SumInt' label 'filterCorrected.fig']);
print(['powerSeries-SumInt' label 'filterCorrected.png'],'-dpng','-r300');
clf();

% plot(power,SumInt,'o');
% xlabel('Excitation Power (mW)');
% %ylabel('overall integrated Intensity (a.u.)');
% ylabel('overall integrated Intensity (a.u.)');
% graphicsSettings;
% savefig('powerSeries-SumInt-filterCorrected-lin.fig');
% print('powerSeries-SumInt-filterCorrected-lin.png','-dpng','-r300');
% clf();

end