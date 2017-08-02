%% How to analyse a 3-Channel dataset
% Data from a 3-Channel homodyne quantum state measurement has its own
% functions for evaluation which will be presented here.

%% Loading a single data-file
% For each dataset you need to specify two files: The file _filenameSIG_
% where you recorded an actual signal, and the file _filenameLO_ where you
% recorded the vacuum. The vacuum signal is used to normalize the computed
% quadratures. Otherwise you would not be able to calculate photon numbers
% later on. X1, X2 and X3 are the calculated quadratures for the
% corresponding channels with their offset already subtracted. It is
% assumed, that X3 has a very slow or no phase-modulation, X1 has a
% modulation of 2w and X2 of w (e.g. w=25Hz).
% piezoSign is 1 or -1 depending on the direction of the first segment of
% the channel with the highest modulation frequency.
[X1,X2,X3, piezoSign] = prepare3ChData(filenameLO, filenameSIG);

%% Plotting Cross-Correlations
% For a first overview of your dataset it can be helpful to plot the 3
% possible cross-correlations between the channels. The following
% cross-correlations originate from a thermal light measurement. Therefore,
% the cross-correlations are the only way to see whether there is a defined
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

%% Compute phase between X1 and X3
% When analyzing a 3 Channel dataset, the final result will be a dataset
% (X,theta), from which it is possible to reconstruct the full Wigner
% function. Therefore, two channels, e.g. X1 & X2 are only used to
% reconstruct the phase _theta_ and selecting the region of interest. The
% third channel, X3, delivers the values for X, after all post-selections
% are done. To get a meaningful phase, we first need to compute the
% relative phase between X1 and X3, which is done here.
% ys is the smoothed crosscorrelation between X1 and X3.
[theta,ys] = computePhase(X1,X3,piezoSign);

%%
% The variable _piezoSign_ is important when you want to put the data from
% different measurements together. It tells the function in which direction
% the piezo actuator started the movement in the first segment. In
% principle, in two different measurements, the acuator could move in
% opposite directions in the first segment, which would lead to wrong phase
% values.

%% Phase Validation
% It is possible that _computePhase_ doesn't reconstruct the phase for
% some piezo segments correctly. Therefore, the following function can be
% used to manually select only the correctly reconstructed segments. The
% segments are plotted consecutively. Then, with a keyboard button press,
% the segment is selected and with a mouse button press rejected.
[X1,X2,X3,theta] = segmentValidation(X1,X2,X3,ys,theta);

%% Select datapoints where two channels are orthogonal
% To reconstruct the phase between the signal field and our local
% oscillator, we need to work only on data points, where two homodyne
% channels are orthogonal to each other. At these points, the
% cross-correlation of the two channels is on average zero. The following
% function selects such quadruplets (O1,O2,O3,oTheta) from
% (X1,X2,X3,theta), where X1 and X2 are orthogonal. Furthermore, there are
% two kinds of orthogonality: pi/2 and 3*pi/2. The function distinguishes
% them and updates the sign of the O2-values accordingly.
[O1,O2,O3,oTheta] = selectOrthogonal(X1,X2,X3,theta,piezoSign,'Plot','plot');

%%
% The default setting selects a range of 5% (Max-Min) around zero. However,
% you can choose your own range by specifying the 'Width' parameter:
%
%   selectOrthogonal(...,'Width',width)

%% 2D Histogram
% It is sometimes useful to calculate and visualize a two dimensional
% histogram of quadrature data. For example of _O1_ and _O2_ from the
% previous function outputs. The amplitudes of these orthogonal channels
% give a point in a phasor diagram and therefore a phase.
[H, binsO1, binsO2] = histogram2D(O1,O2,'plot');

%% Visualizing the selection process with a Movie
% Next, we select from the orthogonally pre-selected dataset (O1,O2,O3)
% only such points, that lie in a specific region of the previously plotted
% histogram. If we plot a histogram of such data for the third channel, O3,
% this should result in a phase-averaged coherent state, when the original
% measurement was performed on a thermal state. The following function does
% this selection and plotting for many different regions and saves the
% results to a Movie. The input variable _theta_ is empty for now.
plot3ChMovie(O1,O2,O3,oTheta,'nomovie');

%%
% The optional parameter _'nomovie'_ prevents the creation of the movie. The
% blue region in the inset illustrates the selected region. And the main
% plot shows the corresponding histogram of the selected O3-values.

%% Selecting an arbitrary region of points
% The following function selects a specific region in a given orthogonal
% 3-Channel dataset. X consists of the selected O3 values and theta of the
% corresponding phase values (already adjusted with the atan2(O2,O1).
[X,theta] = selectRegion(O1,O2,O3,oTheta,'Plot','show');

%% Visualization of density matrix reconstruction
% To get an impression how well the reconstruction algorithm works and how
% fast it converges, it may be useful to visualize this process with a
% movie.plot
[rho100,history] = computeDensityMatrix(X,theta,'Iterations',100);
history = num2cell(history,[1 2]);
plotMovie(@plotRho,history,'Delays',5,'ZLim',[-0.2 0.2]);

%% Helper functions
% If the data in A only exhibits a certain number of different values (e.g.
% 8-bit, ...), you can get the minimum number _nbins_ of bins _bins_.
%
%   [bins, nBins] = minBins(A)
%
% The generic function _plotMovie_ takes any suitable plotting function as
% input to generate the plots. One of them creates a 3D-bar chart of _rho_:
%
%   plotRho(rho);