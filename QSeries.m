function [fitFunctions] = QSeries(O1,O2,O3,oTheta,varargin)
%QSERIES - 
% Input arguments: Selected datapoints from a 3-Channel-Measurement where
% the two channels O1 and O2 are orthogonal. (See in demo3ChannelAnalysis
% how this is done.)
%
% Output arguments: fitFunctions: a struct containing the fit fuctions for
% the Q(postselect) - dependence of <Q>max, <Q^2>max and the uncertainty product. 
%
% Depending on 'Type', a region like a line or a fullcircle with a series of 
% positions Q is postselected from the orthogonal channels phase space,
% and the expectation values of the target channel are computed  and
% plotted.
% You can choose with 'PlotExpectations' whether to plot the expectation
% values for every Q. 
% 'MeasNumber' is the number of the measurement from which the data comes.
% It will be displayed in the plot titles and the filenames.

%% Validate and parse input arguments
p = inputParser;
defaultType = 'Qline';
defaultPlotOpt = 'plot';
defaultMeasNumber = 11;
addParameter(p,'Type',defaultType,@isstr);
addParameter(p,'PlotExpectations',defaultPlotOpt,@isstr);
addParameter(p,'MeasNumber',defaultMeasNumber,@isvector);
parse(p,varargin{:});
c = struct2cell(p.Results);
[measNumber, plotOpt, type] = c{:};

%% Compute expectation values for every Q and plot them

Q = linspace(0,6,25);
width = 0.5;
[ expQmax,  expQ2max, expQ2min, delQmean, Unc, N]=deal(zeros(length(Q),1));

for i = 1:length(Q)
    [XQ,thetaQ] = selectRegion(O1,O2,O3,oTheta,'Type',type,'Position',...
        [Q(i) width],'Plot','hide');        
    [XQdis,thetaQdis]=discretizeTheta(XQ,thetaQ,180);
    [ expQ, ~, expQ2, ~, delQ, ~, meanUnc,~ , meanN, ~] =...
        computeExpectations2( XQdis, thetaQdis,[num2str(measNumber) '-' type '-Q-' ...
        strrep(num2str(Q(i)),'.',',') '-w-' strrep(num2str(width),'.',',')],...
        'Plot', plotOpt );
    expQmax(i) = max(expQ);
    expQ2max(i) = max(expQ2);
    expQ2min(i) = min(expQ2);
    delQmean(i) = mean(delQ);
    Unc(i) = meanUnc;
    N(i) = meanN;
end

%% Fit QPS- dependence of ...

% Fit <Q>max to linear function
fQ = fit(Q', expQmax, 'poly1');

% Fit <Q^2>max to parabolic function
startpoint = [expQ2max(1) (expQ2max(end)-expQ2max(1))/Q(end)^2];
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0],...
               'Upper',[Inf,Inf],...
               'StartPoint',startpoint);
ft = fittype('m*x^2 + b','options',fo);
[fQ2,~] = fit(Q',expQ2max,ft);

% Fit <delQ*delP> to linear function
fUnc = fit(Q', Unc, 'poly1');

fitFunctions = struct('fQ', fQ, 'fQ2', fQ2, 'fUnc', fUnc);

%% Plot expectation values over Q
close all;
hold off;
plot(Q, expQmax, Q, expQ2max, Q, expQ2min, Q, delQmean, Q, Unc, Q, N, 'linewidth', 2);
hold on;
plot(fQ, 'k--');
plot(fQ2, 'k-');
plot(fUnc, 'k-.');

xlabel('Q_{PS}');
%axis([0 max(Q) 0  ]);
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultTextInterpreter','latex');
legend('$<Q>_{max}$', '$<Q^{2}>_{max}$', '$<Q^{2}>_{min}$', '$<\Delta Q>$',...
     '$<\Delta Q \cdot \Delta P>$ ', '$<n>$ ',...
     ['Fit $<Q>_{max} =$' num2str(fQ.p1) '$*Q_{PS} + $' num2str(fQ.p2)],...
     ['Fit $<Q^2>_{max} =$' num2str(fQ2.m) '$*Q_{PS}^2 + $' num2str(fQ2.b)],...
     ['Fit $<\Delta Q \cdot \Delta P>=$' num2str(fUnc.p1) '$*Q_{PS} + $' num2str(fUnc.p2)],...
     'location', 'best');
title([num2str(measNumber) '-' type '-Width-' num2str(width)]);

    % Save plot
print([num2str(measNumber) '-' type '-w-' strrep(num2str(width),'.',',') ...
    '-Qdependence'], '-dpng');

hold off;

end