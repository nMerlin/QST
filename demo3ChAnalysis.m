%% How to analyse a 3-Channel dataset
% Data from a 3-Channel homodyne quantum state measurement has its own
% functions for evaluation which will be presented here.

%% Loading a single data-file
% For each dataset you need to specify two files: The file _filenameSIG_
% where you recorded an actual signal, and the file _filenameLO_ where you
% recorded the vacuum. The vacuum signal is used to normalize the computed
% quadratures. Otherwise you would not be able to calculate photon numbers
% later on. X1, X2 and X3 are the calculated quadratures for the
% corresponding channels with their offset already subtracted.
[X1,X2,X3] = prepare3ChData(filenameLO, filenameSIG);

%% Plotting Cross-Correlations
% For a first overview of your dataset it can be helpful to plot the 3
% possible cross-correlations between the channels. The following
% cross-correlations originate from a thermal light measurement. Therefore,
% the cross-correlations are the only way to see wether there is a defined
% phase between your channels or not. Per default, the plot is limited to
% two consecutive piezo segments.
plotCrossCorrelation(X1,X2,X3);

%%
% The previous plotting function makes use of the function
%
%   ys = smoothCrossCorr(Xa,Xb,varargin)
%
% which generates a smoothed version of the vector
%
%   Xa.*Xb
%
% with the help of cubic smoothing splines.

%% Select datapoints where two channels are orthogonal
% To reconstruct the phase between the signal field and our local
% oscillator, we need to work only on data points, where two homodyne
% channels are orthogonal to each other. At these points, the
% cross-correlation of the two channels is on average zero. The following
% function selects such triples (O1,O2,O3) from (X1,X2,X3), where X1 and X2
% are orthogonal.
[O1,O2,O3] = selectOrthogonal(X1,X2,X3,'plot');

%%
% The default setting selects a range of 5% (Max-Min) around zero. However,
% you can choose your own range by specifying the _ORTH_WIDTH_ parameter:
%
%   selectOrthogonal(X1,X2,X3,'noplot',ORTH_WIDTH)