function [] = comparisonPolarizationsFromNPhotons()

% get data
cd('01-l4-10-l2-0-GT-0');
load('AverageNandG2-X-time-weight-yes.mat');
power1 = Is;
n1 = nAvs;
g2_1 = g2Avs;
g2Stds_1 = g2Stds;
cd('..');
cd('02-l4-10-l2-22.5-GT-0');
load('AverageNandG2-X-time-weight-yes.mat');
power2 = Is;
n2 = nAvs;
g2_2 = g2Avs;
g2Stds_2 = g2Stds;
cd('..');
cd('03-l4-10-l2-45-GT-0');
load('AverageNandG2-X-time-weight-yes.mat');
power3 = Is;
n3 = nAvs;
g2_3 = g2Avs;
g2Stds_3 = g2Stds;
cd('..');
cd('04-l4-10-l2-67.5-GT-0');
load('AverageNandG2-X-time-weight-yes.mat');
power4 = Is;
n4 = nAvs;
g2_4 = g2Avs;
g2Stds_4 = g2Stds;
cd('..');
cd('05-l4-55-l2-0-GT-0');
load('AverageNandG2-X-time-weight-yes.mat');
power5 = Is;
n5 = nAvs;
g2_5 = g2Avs;
g2Stds_5 = g2Stds;
cd('..');
cd('06-l4-325-l2-0-GT-0');
load('AverageNandG2-X-time-weight-yes.mat');
power6 = Is;
n6 = nAvs;
g2_6 = g2Avs;
g2Stds_6 = g2Stds;
cd('..');

%% computation of degree of polarisation
P1 = (n1 - n3)./(n1 + n3); 
P2 = (n2 - n4)./(n2 + n4); 
P3 = (n5 - n6)./(n5 + n6); 
Ptot = sqrt(P1.^2 + P2.^2 + P3.^2);
n_tot = n1+n2+n3+n4+n5+n6;

%% plot of all n
loglog(power1,n1,'o-','Linewidth',2); hold on;
loglog(power2,n2,'o-','Linewidth',2);
loglog(power3,n3,'o-','Linewidth',2);
loglog(power4,n4,'o-','Linewidth',2);
loglog(power5,n5,'o-','Linewidth',2);
loglog(power6,n6,'o-','Linewidth',2);
legend('l/4 10°, l/2 0°', 'l/4 10°, l/2 22.5°','l/4 10°, l/2 45°','l/4 10°, l/2 67.5°',...
    'l/4 55°, l/2 0°','l/4 325°, l/2 0°','location','northwest');
xlabel('Excitation Power (mW)');
ylabel('Photon Number');
graphicsSettings;
xlim([1 100]);
savefig('comparison-n-all.fig');
print('comparison-n-all.png','-dpng','-r300');
clf();


%% plot of pol. degree 
semilogx(power1,P1,'o-','Linewidth',2); hold on;
semilogx(power2,P2,'o-','Linewidth',2);
semilogx(power3,P3,'o-','Linewidth',2);
semilogx(power4,Ptot,'o-','Linewidth',2);
legend('P1: l/2 0° vs 45°', 'P2: l/2 22.5° vs 67.5°',...
    'P3: l/4 55° vs 325°', 'P_{tot}','location','northwest');
xlabel('Excitation Power (mW)');
ylabel('Polarisation Degree');
ylim([-1 1]);
graphicsSettings;
xlim([1 100]);
savefig('comparison-polDegree-n-all.fig');
print('comparison-polDegree-n-all.png','-dpng','-r300');
clf();


%% plot of sum n
loglog(power1,n_tot,'o','Linewidth',2);
legend('Sum of all polarizations','location','northwest');
xlabel('Excitation Power (mW)');
ylabel('Photon Number');
graphicsSettings;
xlim([1 100]);
savefig('comparison-n-all-sum.fig');
print('comparison-n-all-sum.png','-dpng','-r300');
clf();

%% plot of all g2
hAx=axes;
hAx.XScale='log';
%xlim([minDecade maxDecade])
hold all;
% errorbar(power1,g2_1,g2Stds_1,'o-','Linewidth',2);
% errorbar(power2,g2_2,g2Stds_2,'o-','Linewidth',2);
% errorbar(power3,g2_3,g2Stds_3,'o-','Linewidth',2);
% errorbar(power4,g2_4,g2Stds_4,'o-','Linewidth',2);
% errorbar(power5,g2_5,g2Stds_5,'o-','Linewidth',2);
% errorbar(power6,g2_6,g2Stds_6,'o-','Linewidth',2);
plot(power1,g2_1,'o-','Linewidth',2);
plot(power2,g2_2,'o-','Linewidth',2);
plot(power3,g2_3,'o-','Linewidth',2);
plot(power4,g2_4,'o-','Linewidth',2);
plot(power5,g2_5,'o-','Linewidth',2);
plot(power6,g2_6,'o-','Linewidth',2);
legend('l/4 10°, l/2 0°', 'l/4 10°, l/2 22.5°','l/4 10°, l/2 45°','l/4 10°, l/2 67.5°',...
    'l/4 55°, l/2 0°','l/4 325°, l/2 0°','location','northwest');
xlabel('Excitation Power (mW)');
ylabel('$ \langle g^{(2)}(0) \rangle $','Interpreter','latex');
ylim([0.5 2.5]);
xlim([1 100]);
graphicsSettings;
savefig('comparison-g2-all.fig');
print('comparison-g2-all.png','-dpng','-r300');
clf();

end