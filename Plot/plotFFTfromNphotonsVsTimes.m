function [] = plotFFTfromNphotonsVsTimes(ada,times,filename,varargin)
%% PLOTFFT plots the fourier transformation of the time dependent photon number ada.
%
%   Input Arguments:
%       ada: time dependent photon numbers, that are computed with 
%       [g2vec, ada, times] = g2(X,nResolution, varargin) from the quadratures X
%       times: time vector
%       filename: name for the saved figure

%% Validate and parse input arguments
p = inputParser;
defaultAxis = [0 7000 0 1];
addParameter(p,'Axis',defaultAxis);
parse(p,varargin{:});
c = struct2cell(p.Results);
[Ax] = c{:};

    T = mean(diff(times));
    samplerate = 1/T;

    % Calculating
    L = length(ada);
    Y = fft(ada);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = samplerate*(0:(L/2))/L;

    % Plotting
    [pks,locs] = findpeaks(P1,f,'MinPeakProminence',max(P1)/20);
    findpeaks(P1,f,'MinPeakProminence',max(P1)/20);
    text(locs+10,pks,num2str(locs'));
    title('Single-Sided Amplitude Spectrum');
    xlabel('f (Hz)');
    ylabel('|P1(f)|');
    axis(Ax);

    % Saving 
    savefig([filename '-FFT.fig']);
    print([filename '-FFT.png'],'-dpng','-r300');
end