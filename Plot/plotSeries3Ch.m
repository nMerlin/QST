function plotSeries3Ch(varargin)
%PLOTSERIES3CH Plot results from SERIES3CH time series

delays = T.Delay;
varX = T.meanVarX;
plot(delays,varX,'o');
fo = fitoptions('Method','NonlinearLeastSquares', ...
    'StartPoint',[-5 0 600 8.3]);
ft = fittype('a*exp(-(x-b)^2/t^2)+c','options',fo);
res = fit(delays,varX,ft);
hold on;
plot(-2000:1:2000,res(-2000:1:2000),'r');
xlim([-2000 2000]);
xlabel('Delay (fs)');
ylabel('Postselected Variance');

% Lower Bound
minVar = compute3ChLimit(T.nX1,T.nX2,T.nX3);
plot(delays,minVar,'x');

% Upper Bound
plot(delays,T.nX3+0.5,'*');

legend(['Minimum: ',num2str(min(varX))],'Fit', ...
    ['Limit ca. ',num2str(max(minVar))],'Upper Bound', ...
    'Location','SouthEast');

hold off;

end

