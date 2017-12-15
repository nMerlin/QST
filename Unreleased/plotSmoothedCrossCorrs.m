function plotSmoothedCrossCorrs(Xa, Xb, params,varargin)
%PLOTSMOOTHEDCROSSCORRS Plot of smoothed cross-correlation for different
%parameters.
%
% Input Arguments:
%   Xa - matrix of quadratures in first channel
%   Xb - matrix of quadratures in second channel
%   params - vector of all smoothing parameters to plot

%% Constants
x = (1:length(Xa(:)))*1/75.4e6*1e3;

%% Validate and parse input arguments
p = inputParser;
defaultRange = [min(x) max(x)];
addParameter(p,'Range',defaultRange,@isvector);
parse(p,varargin{:});
c = struct2cell(p.Results);
[range] = c{:};

nParams = length(params);
for k=1:nParams
    ys = smoothCrossCorr(Xa,Xb,'Param',params(k));
    if k~=1
        hold on;
    end
    plot(x,ys);
end
set(gca,'XLim',range);
legend('1e-10','1e-11','1e-12','1e-13','1e-14','1e-15');
title('Simulated Data: Smoothed Cross-Correlations');
xlabel('Time [ms]');
ylabel('Smoothing Splines of X1*X2');

end

