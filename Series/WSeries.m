function [] = WSeries(O1,O2,O3,oTheta,varargin)
%WSERIES - 
% Input arguments: Selected datapoints from a 3-Channel-Measurement where
% the two channels O1 and O2 are orthogonal. (See in demo3ChannelAnalysis
% how this is done.)
%
% Output arguments: fitFunctions: a struct containing the fit fuctions for
% the w(postselect) - dependence of <Q>max, <Q^2>max and the uncertainty product. 
%
% Depending on 'Type', a region like a line or a fullcircle with a series of 
% widths w is postselected from the orthogonal channels phase space,
% and the expectation values of the target channel are computed  and
% plotted.
% You can choose with 'PlotExpectations' whether to plot the expectation
% values for every w. 
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

w = linspace(0.01,5,30);
Q = 2.5;

[ expQmax,  expQ2max, expQ2min, delQmean, Unc, N, nEvents]=deal(zeros(length(w),1));

for i = 1:length(w)
    [XQ,thetaQ] = selectRegion(O1,O2,O3,oTheta,'Type',type,'Position',...
        [Q w(i)],'Plot','hide');        
    XQ = XQ(~isnan(thetaQ));
    thetaQ=thetaQ(~isnan(thetaQ)); %remove X and theta values where theta is nan.
    [XQdis,thetaQdis]=discretizeTheta(XQ,thetaQ,180);
    [ expQ, ~, expQ2, ~, delQ, ~, meanUnc,~ , meanN, ~] =...
        computeExpectations2( XQdis, thetaQdis,[num2str(measNumber) '-' type '-Q-' ...
        strrep(num2str(Q),'.',',') '-w-' strrep(num2str(w(i)),'.',',')],...
        'Plot', plotOpt );
    expQmax(i) = max(expQ);
    expQ2max(i) = max(expQ2);
    expQ2min(i) = min(expQ2);
    delQmean(i) = mean(delQ);
    Unc(i) = meanUnc;
    N(i) = meanN;
    nEvents(i) = length(XQ);
end

nEvents = nEvents / max(nEvents); %Normalize nEvents 

%% Plot expectation values over w
close all;
hold off;
plot(w, expQmax, w, expQ2max, w, expQ2min, w, delQmean, w, Unc, w, N, w, nEvents,'linewidth', 2);

xlabel('w_{PS}');
%axis([0 max(Q) 0  ]);
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultTextInterpreter','latex');
legend('$<Q>_{max}$', '$<Q^{2}>_{max}$', '$<Q^{2}>_{min}$', '$<\Delta Q>$',...
     '$<\Delta Q \cdot \Delta P>$ ', '$<n>$ ', 'norm. N of Events', 'location', 'best');
title([num2str(measNumber) '-' type '-$Q_{PS}$-' num2str(Q)]);

    % Save plot
print([num2str(measNumber) '-' type '-Q_{PS}-' strrep(num2str(Q),'.',',') ...
    '-wDependence'], '-dpng');

end