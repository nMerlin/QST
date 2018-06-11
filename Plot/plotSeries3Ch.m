function plotSeries3Ch(T,varargin)
%PLOTSERIES3CH Plot results from SERIES3CH table T.

%% Validate and parse input arguments
p = inputParser;
defaultSave = false;
addParameter(p,'Save',defaultSave,@islogical);
parse(p,varargin{:});
c = struct2cell(p.Results);
[saveflag] = c{:};

%% Create plot
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

% Lower Bound
minVar = compute3ChLimit(T.nX2,T.nX3,T.nX1);
plot(delays,minVar,'x');

% Upper Bound
plot(delays,T.nX1+0.5,'*');

% Labels
xlabel('Delay (fs)');
ylabel('Postselected Variance');
legend(['Minimum: ',num2str(min(varX))],'Fit', ...
    ['Limit ca. ',num2str(max(minVar))],'Upper Bound', ...
    'Location','SouthEast');
title('3-Channel Series');
hold off;

%% Write figure to file
if saveflag
    saveA5Landscape([datestr(date,'yyyy-mm-dd-'),'plotSeries3Ch']);
end

end

