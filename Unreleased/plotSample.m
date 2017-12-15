function [] = plotSample(sample,ys)
%PLOTSAMPLE Plot sample data.

xSample = (1:length(sample))*1/(75.4e6)*1e3;

plot(xSample,sample,'.');
hold on;
h = plot(xSample,ys,'r-');
xlabel('Time [ms]');
ylabel('X1*X2/Smoothing Spline');
title('Sample B: Smoothed Cross-Correlations');
set(gca,'XLim',[min(xSample) max(xSample)]);
set(gca,'YLim',[-30 30]);
set(h,'LineWidth',5);
legend('X1*X2','Smoothing Spline 1e-15');
hold off;

end

