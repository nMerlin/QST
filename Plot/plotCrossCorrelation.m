function plotCrossCorrelation(X1, X2, X3, varargin)
%PLOTCROSSCORRELATION Plot all three possible crosscorrelations
%
%   (X1,X2,X3) is a 3-Channel dataset, prepared by PREPARE3CHDATA (i.e.
%   reshaped into piezo-segments with offset already removed)
%
%   Optional input arguments:
%   plotCrossCorrelation(X1,X2,X3,N_SEGMENTS): plot N_SEGMENTS piezo
%       segments, default is 2

%% Handle optional input arguments and default values
nVarargin = length(varargin);
optArgs = {2};
optArgs(1:nVarargin) = varargin;
[N_SEGMENTS] = optArgs{:};

%% Compute smoothed cross-correlations
ys12 = smoothCrossCorr(X1(:,:,1:N_SEGMENTS),X2(:,:,1:N_SEGMENTS));
ys13 = smoothCrossCorr(X1(:,:,1:N_SEGMENTS),X3(:,:,1:N_SEGMENTS));
ys23 = smoothCrossCorr(X2(:,:,1:N_SEGMENTS),X3(:,:,1:N_SEGMENTS));

%% Plot
%[nPoints,~] = size(ys12);
hold on;
plot(ys12(:),'linewidth',3);
plot(ys13(:),'linewidth',3);
plot(ys23(:),'linewidth',3);
hold off;
%set(gca,'XLim',[min(x) max(x)]);
title('Smoothed Cross-Correlations');
legend('X1*X2','X1*X3','X2*X3');

end

