function [ noise, nResolution ] = compareG2s( X, filename )
%COMPAREG2S Compares the noise of g2 for several integration times
%   Detailed explanation goes here

samplingrate = 75.4; % MHz

%%% Use Smoothingspline-Fit to get a "signal"-trace
vecG2 = g2(X,2000);
xAxis = (1:length(vecG2))'*2000*1/samplingrate;
signalFit = fit(xAxis,vecG2,'smoothingspline', ...
    'smoothingparam',0.0000005);
clf;
plot(xAxis,vecG2,'.');
hold on;
plot(signalFit);
set(gca,'YLim',[0.9 2.8],'XLim',[min(xAxis) max(xAxis)]);
xlabel('Time [{\mu}s]');
ylabel('g^{(2)}(0)');
legend('26.5 {\mu}s per Point: g^{(2)}_{meas}', ...
    'Smoothingspline: g^{(2)}_{spline}');
title('Noise vs. Time Constant for g^{(2)}');

%%% Compute the standard deviation for each value in "nResolution"
nResolution = 100:100:2000;
noise = zeros(length(nResolution),1);
for k = 1:length(nResolution)
    vecG2 = g2(X,nResolution(k));
    xAxis = (1:length(vecG2))*nResolution(k)/samplingrate;
    signal = signalFit(xAxis);
    noise(k) = sqrt(var(vecG2-signal));
end

%%% Inset with noise behavior
insetAx = axes('Parent',gcf,'Position',[0.23 0.6 0.3 0.25]);
noiseX = nResolution*1/samplingrate;
loglog(noiseX,noise,'o');
hold on;
title('$$\sigma = \sqrt{VAR({g^{(2)}_{meas} - g^{(2)}_{spline}})}$$', ...
    'Interpreter','latex');
set(insetAx,'FontSize',8,'XLim',[min(noiseX) max(noiseX)], ...
    'YLim',[min(noise) max(noise)],'XTick',[1 2 5 10 20]);
ylabel('\sigma (Std. Dev.)');
xlabel('Time {\Delta}t per point [{\mu}s]');

%%% Linear fit
p = polyfit(log(noiseX)',log(noise),1);
loglog(noiseX,exp(p(2))*noiseX.^p(1));
fitString = ['{\sigma} = ',num2str(exp(p(2)),2),'{\cdot}','{\Delta}t^{',...
    num2str(p(1),2),'}'];
textX = exp((0.58352*log(max(noiseX))+log(min(noiseX)))/(1+0.58352));
textY = exp((4.882475*log(max(noise))+log(min(noise)))/(1+4.882475));
text(textX,textY,fitString);
hold off;

%%% Save figure
print(strcat('G2-',filename,'.png'), '-dpng');
end

