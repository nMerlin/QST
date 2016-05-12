function [  ] = plotWaveformsLecroy( waveformdata )
% Plot waveforms from LeCroy oscilloscope on top of each other
s = size(waveformdata);
x = 0:200/(length(waveformdata)-1):200;

for i=1:s(2)
    y = waveformdata(:,i);
    y = y - mean(y);
    plot(x,y,'b');
    hold on;
end

xlabel('time [ns]');
ylabel('voltage [V]');
print('plottedWaveforms.png','-r300','-dpng');

end