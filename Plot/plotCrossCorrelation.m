function [A12,A13,A23] = plotCrossCorrelation(X1, X2, X3, filename, varargin)
%PLOTCROSSCORRELATION Plot all three possible crosscorrelations
%
%   (X1,X2,X3) is a 3-Channel dataset, prepared by PREPARE3CHDATA (i.e.
%   reshaped into piezo-segments with offset already removed)
%
%   Optional input arguments:
%   plotCrossCorrelation(X1,X2,X3,N_SEGMENTS): plot N_SEGMENTS piezo
%       segments, default is 2

%% Validate and parse input arguments
p = inputParser;
defaultNSegments = 2;
addParameter(p,'NSegments',defaultNSegments,@isnumeric);
parse(p,varargin{:});
c = struct2cell(p.Results);
[nsegments] = c{:};

%% Compute smoothed cross-correlations
ys12 = smoothCrossCorr(X1(:,:,1:nsegments),X2(:,:,1:nsegments));
ys13 = smoothCrossCorr(X1(:,:,1:nsegments),X3(:,:,1:nsegments));
ys23 = smoothCrossCorr(X2(:,:,1:nsegments),X3(:,:,1:nsegments));

    function [A] = vis(y)  %compute the visibility
        y = y - min(min(y)); % The function is symmetric around zero and we want it to become positive
        [maxpks, ~] = findpeaks(y(:),'minPeakProminence',0.2*abs(max(max(y))));
        [minpks, ~] = findpeaks(-y(:),'minPeakProminence',0.2*abs(max(max(-y))));
        minpks = -minpks;
        maxP = mean(maxpks);
        minP = mean(minpks);
        %A = (maxP - minP)/(maxP + minP);
        A = (maxP - minP);
    end


A12 = vis(ys12);
A13 = vis(ys13);
A23 = vis(ys23);

%% Plot
plot(ys12(:),'linewidth',3);
hold on;
plot(ys13(:),'linewidth',3);
plot(ys23(:),'linewidth',3);
hold off;
set(gca,'XLim',[1 length(ys12(:))]);
title('Smoothed Cross-Correlations');
legend('X1*X2','X1*X3','X2*X3');
graphicsSettings;
savefig([filename '-Cross-Correlations' '.fig']);
print([filename '-Cross-Correlations' '.png'],'-dpng','-r300');
clf();

end

