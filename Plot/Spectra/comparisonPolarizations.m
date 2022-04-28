function [] = comparisonPolarizations()

% get data
% M = xlsread('powerSeries-filterCorrected.xls','A2:L50');
M = xlsread('powerSeries-l4-10-l2-5-filterCorrected.xls','A2:L50');
power1 = M(:,8);
modeInt1 = M(:,9);
SumInt1 = M(:,10);
%remove NaNvalues
power1 = power1(~isnan(power1));
modeInt1 = modeInt1(~isnan(modeInt1));
SumInt1 = SumInt1(~isnan(SumInt1));
M = xlsread('powerSeries-l4-10-l2-27.5-filterCorrected.xls','A2:L50');
power2 = M(:,8);
modeInt2 = M(:,9);
SumInt2 = M(:,10);
%remove NaNvalues
power2 = power2(~isnan(power2));
modeInt2 = modeInt2(~isnan(modeInt2));
SumInt2 = SumInt2(~isnan(SumInt2));
M = xlsread('powerSeries-l4-10-l2-50-filterCorrected.xls','A2:L50');
power3 = M(:,8);
modeInt3 = M(:,9);
SumInt3 = M(:,10);
%remove NaNvalues
power3 = power3(~isnan(power3));
modeInt3 = modeInt3(~isnan(modeInt3));
SumInt3 = SumInt3(~isnan(SumInt3));
M = xlsread('powerSeries-l4-10-l2-72.5-filterCorrected.xls','A2:L50');
power4 = M(:,8);
modeInt4 = M(:,9);
SumInt4 = M(:,10);
%remove NaNvalues
power4 = power4(~isnan(power4));
modeInt4 = modeInt4(~isnan(modeInt4));
SumInt4 = SumInt4(~isnan(SumInt4));
M = xlsread('powerSeries-l4-55-l2-5-filterCorrected.xls','A2:L50');
power5 = M(:,8);
modeInt5 = M(:,9);
SumInt5 = M(:,10);
%remove NaNvalues
power5 = power5(~isnan(power5));
modeInt5 = modeInt5(~isnan(modeInt5));
SumInt5 = SumInt5(~isnan(SumInt5));
M = xlsread('powerSeries-l4-325-l2-5-filterCorrected.xls','A2:L50');
power6 = M(:,8);
modeInt6 = M(:,9);
SumInt6 = M(:,10);
%remove NaNvalues
power6 = power6(~isnan(power6));
modeInt6 = modeInt6(~isnan(modeInt6));
SumInt6 = SumInt6(~isnan(SumInt6));

%% computation of degree of polarisation
modeInt1 = modeInt1(1:length(power2));
SumInt1 = SumInt1(1:length(power2));
power1 = power1(1:length(power2));
P1_Mode = (modeInt1 - modeInt3)./(modeInt1 + modeInt3); 
P2_Mode = (modeInt2 - modeInt4)./(modeInt2 + modeInt4); 
P3_Mode = (modeInt5 - modeInt6)./(modeInt5 + modeInt6); 
Ptot_Mode = sqrt(P1_Mode.^2 + P2_Mode.^2 + P3_Mode.^2);
P1_Sum = (SumInt1 - SumInt3)./(SumInt1 + SumInt3); 
P2_Sum = (SumInt2 - SumInt4)./(SumInt2 + SumInt4); 
P3_Sum = (SumInt5 - SumInt6)./(SumInt5 + SumInt6); 
Ptot_Sum = sqrt(P1_Sum.^2 + P2_Sum.^2 + P3_Sum.^2);
SumInt_tot = SumInt1+SumInt2+SumInt3+SumInt4+SumInt5+SumInt6;
modeInt_tot = modeInt1+modeInt2+modeInt3+modeInt4+modeInt5+modeInt6;

%% plot of all modeInt
loglog(power1,modeInt1,'o-','Linewidth',2); hold on;
loglog(power2,modeInt2,'o-','Linewidth',2);
loglog(power3,modeInt3,'o-','Linewidth',2);
loglog(power4,modeInt4,'o-','Linewidth',2);
loglog(power5,modeInt5,'o-','Linewidth',2);
loglog(power6,modeInt6,'o-','Linewidth',2);
legend('l/4 10°, l/2 5°', 'l/4 10°, l/2 27.5°','l/4 10°, l/2 50°','l/4 10°, l/2 72.5°',...
    'l/4 55°, l/2 5°','l/4 325°, l/2 5°','location','northwest');
xlabel('Excitation Power (mW)');
ylabel('Mode Intensity (counts/time)');
graphicsSettings;
savefig('comparison-modeInt-all-filterCorrected.fig');
print('comparison-modeInt-all-filterCorrected.png','-dpng','-r300');
clf();

%% plot of all Integrated
loglog(power1,SumInt1,'o-','Linewidth',2,'MarkerFaceColor','auto'); hold on;
loglog(power2,SumInt2,'o-','Linewidth',2,'MarkerFaceColor','auto');
loglog(power3,SumInt3,'o-','Linewidth',2,'MarkerFaceColor','auto');
loglog(power4,SumInt4,'o-','Linewidth',2,'MarkerFaceColor','auto');
loglog(power5,SumInt5,'o-','Linewidth',2,'MarkerFaceColor','auto');
loglog(power6,SumInt6,'o-','Linewidth',2,'MarkerFaceColor','auto');
legend('l/4 10°, l/2 5°', 'l/4 10°, l/2 27.5°','l/4 10°, l/2 50°','l/4 10°, l/2 72.5°',...
    'l/4 55°, l/2 5°','l/4 325°, l/2 5°','location','northwest');
xlabel('Excitation Power (mW)');
ylabel('overall integrated Intensity (counts/time)');
graphicsSettings;
savefig('comparison-SumInt-all-filterCorrected.fig');
print('comparison-SumInt-all-filterCorrected.png','-dpng','-r300');
clf();

%% plot of pol. degree mode
semilogx(power1,P1_Mode,'o-','Linewidth',2); hold on;
semilogx(power2,P2_Mode,'o-','Linewidth',2);
semilogx(power3,P3_Mode,'o-','Linewidth',2);
semilogx(power4,Ptot_Mode,'o-','Linewidth',2);
legend('P1: l/2 5° vs 50°', 'P2: l/2 27.5° vs 72.5°',...
    'P3: l/4 55° vs 325°', 'P_{tot}','location','northwest');
xlabel('Excitation Power (mW)');
ylabel('Polarisation Degree of mode');
ylim([-1 1]);
graphicsSettings;
savefig('comparison-polDegree-modeInt-all-filterCorrected.fig');
print('comparison-polDegree-modeInt-all-filterCorrected.png','-dpng','-r300');
clf();

%% plot of pol. degree integrated
semilogx(power1,P1_Sum,'o-','Linewidth',2); hold on;
semilogx(power2,P2_Sum,'o-','Linewidth',2);
semilogx(power3,P3_Sum,'o-','Linewidth',2);
semilogx(power4,Ptot_Sum,'o-','Linewidth',2);
legend('P1: l/2 5° vs 50°', 'P2: l/2 27.5° vs 72.5°',...
    'P3: l/4 55° vs 325°', 'P_{tot}','location','northwest');
xlabel('Excitation Power (mW)');
ylabel('Polarisation Degree overall integrated');
ylim([-1 1]);
graphicsSettings;
savefig('comparison-polDegree-SumInt-all-filterCorrected.fig');
print('comparison-polDegree-SumInt-all-filterCorrected.png','-dpng','-r300');
clf();

%% plot of sum modeInt
loglog(power1,modeInt_tot,'o','Linewidth',2);
legend('Sum of all polarizations','location','northwest');
xlabel('Excitation Power (mW)');
ylabel('Mode Intensity (counts/time)');
graphicsSettings;
savefig('comparison-modeInt-all-sum-filterCorrected.fig');
print('comparison-modeInt-all-sum-filterCorrected.png','-dpng','-r300');
clf();

%% plot of sum integrated
loglog(power1,SumInt_tot,'o','Linewidth',2);
legend('Sum of all polarizations','location','northwest');
xlabel('Excitation Power (mW)');
ylabel('overall integrated (counts/time)');
graphicsSettings;
savefig('comparison-SumInt-all-sum-filterCorrected.fig');
print('comparison-SumInt-all-sum-filterCorrected.png','-dpng','-r300');
clf();
end