% get unconverted data
suffix = '-logFit-noLevel';

M = xlsread(['Unconverted-MiraPowerSeries' suffix '.xls'],'B2:C40');
Powers = M(1:end-1,1);
time = M(1:end-1,2);

% get converted data
M = xlsread(['converted-MiraPowerSeries' suffix '.xls'],'B2:C40');
ConvPowers = M(1:end-1,1);
ConvTime = M(1:end-1,2);

% %% plot 
% plot(Powers, time,'o','LineWidth',2);
% title('OPO power - converted light');
% ylabel('decay time (ps)');
% xlabel('OPO Power (mW)');
% graphicsSettings;
% savefig('Converted-OPOSeries.fig');
% print('Converted-OPOSeries.png','-dpng','-r300');


%% plot unconverted
semilogx(Powers, time,'o','LineWidth',2);
title('MIRA power - unconverted light');
ylabel('decay time (ps)');
xlabel('MIRA Power (mW)');
graphicsSettings;
savefig(['Unconverted-MIRASeries' suffix '.fig']);
print(['Unconverted-MIRASeries' suffix '.png'],'-dpng','-r300');

%% plot converted
semilogx(ConvPowers, ConvTime,'o','LineWidth',2);
title('MIRA power - converted light');
ylabel('decay time (ps)');
xlabel('MIRA Power (mW)');
graphicsSettings;
savefig(['Converted-MIRASeries' suffix '.fig']);
print(['Converted-MIRASeries' suffix '.png'],'-dpng','-r300');

% plot factor
semilogx(ConvPowers, ConvTime./time,'o','LineWidth',2);
title('MIRA power - converted time / unconverted time');
ylabel('\tau_{con}/\tau_{un}');
xlabel('MIRA Power (mW)');
graphicsSettings;
savefig(['ConvUnc-MIRASeries' suffix '.fig']);
print(['ConvUnc-MIRASeries' suffix '.png'],'-dpng','-r300');

