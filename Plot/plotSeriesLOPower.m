function plotSeriesLOPower(filename)
%PLOTSERIESLOPOWER Plots LO power series evaluated with seriesLOPower

T = readtable(filename);
errorbar(T.PowerLO(1:21),T.OutputPowerHD(1:21), ...
    T.StandardDeviation(1:21),'o','Color','black','MarkerFaceColor', ...
    'black','MarkerSize',4);
hold on;
errorbar(T.PowerLO(22:41),T.OutputPowerHD(22:41), ...
    T.StandardDeviation(22:41),'o','Color',[.6 .6 .6],'MarkerFaceColor',...
    [.6 .6 .6],'MarkerSize',4);
p = polyfit(T.PowerLO,T.OutputPowerHD,1);
f = polyval(p,T.PowerLO);
plot(T.PowerLO,f,'Color','red','LineWidth',1);
xlim([min(T.PowerLO),max(T.PowerLO)]);
xlabel('LO Power (mW)');
ylabel('BD Power (arb. units)');
legend('First Series','Second Series','Linear Fit','Location','northwest');
set(gca,'FontSize',14);
hold off;

fig = gcf;
set(fig,'Color','w');
fig.Units = 'centimeters';
fig.Position = [1,1,21,7.4];
fig.PaperPositionMode = 'auto';
fig.PaperOrientation = 'landscape';
export_fig 'seriesLOPowerHQ.pdf'

end

