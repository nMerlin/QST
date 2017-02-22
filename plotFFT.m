function [  ] = plotFFT( signal, samplerate )
%PLOTFFT taken from https://de.mathworks.com/help/matlab/ref/fft.html
%   SIGNAL: Recorded Spectraum
%   SAMPLERATE: Sampling Frequency

T = 1/samplerate;
L = length(signal);
Y = fft(signal);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = samplerate*(0:L/2+1);

plot(f,P1);
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

end

