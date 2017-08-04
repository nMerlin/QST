function QSeries(O1,O2,O3,oTheta,varargin)
%QSERIES - 
% Input arguments: Selected datapoints from a 3-Channel-Measurement where
% the two channels O1 and O2 are orthogonal. (See in demo3ChannelAnalysis
% how this is done.)
%
% Depending on 'Type', a region like a line or a fullcircle with a series of 
% positions Q is postselected from the orthogonal channels phase space,
% and the expectation values of the target channel are computed  and
% plotted.

%% Validate and parse input arguments
p = inputParser;
defaultType = 'Qline';
addParameter(p,'Type',defaultType,@isstr);
parse(p,varargin{:});
c = struct2cell(p.Results);
[type] = c{:};

%% Compute expectation values for every Q and plot them

Q = linspace(0,6,25);
width = 0.5;
[ expQmax,  expQ2max, expQ2min, delQmean, Unc, N]=deal(zeros(length(Q),1));

for i = 1:length(Q)
    [XQ,thetaQ] = selectRegion(O1,O2,O3,oTheta,'Type',type,'Position',...
        [Q(i) width],'Plot','hide');        
    [XQdis,thetaQdis]=discretizeTheta(XQ,thetaQ,180);
    [ expQ, ~, expQ2, ~, delQ, ~, meanUnc,~ , meanN, ~] =...
        computeExpectations2( XQdis, thetaQdis,['11-' type '-Q-' ...
        strrep(num2str(Q(i)),'.',',') '-w-' strrep(num2str(width),'.',',')] );
    expQmax(i) = max(expQ);
    expQ2max(i) = max(expQ2);
    expQ2min(i) = min(expQ2);
    delQmean(i) = mean(delQ);
    Unc(i) = meanUnc;
    N(i) = meanN;
end

%% Plot expectation values over Q
close all;
hold off;
plot(Q, expQmax, Q, expQ2max, Q, expQ2min, Q, delQmean, Q, Unc, Q, N, 'linewidth', 2);

xlabel('Q_{PS}');
%axis([0 max(Q) 0  ]);
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultTextInterpreter','latex');
legend('$<Q>_{max}$', '$<Q^{2}>_{max}$', '$<Q^{2}>_{min}$', '$<\Delta Q>$',...
     '$<\Delta Q \cdot \Delta P>$ ', '$<n>$ ', 'location', 'best');
title(['11-' type '-Width-' num2str(width)]);

    % Save plot
print(['11-' type '-w-' strrep(num2str(width),'.',',') '-Qdependence'], '-dpng');

end