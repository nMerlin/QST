function [] = plotNPhotonsAndCrossCorr(X1,X2,piezoSign,filename)

ys = smoothCrossCorr(X1,X2,'Type','spline','Param',1e-14);
 
%% compute photon number
n =  photonNumberVector(X1);

%%Plot
ysR = ys(:);
%plot(n(1:15e5));
plot(n(:));
hold on;
plotysR(:);
%plot(ysR(1:15e5));
legend('n','X1*X2');
graphicsSettings;
print([filename '-nPhotons-CrosscorrelationX1X2.png'],'-dpng','-r300');
savefig([filename '-nPhotons-CrosscorrelationX1X2.fig']);
clf();

end