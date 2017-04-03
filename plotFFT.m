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
f = samplerate*(0:(L/2))/L;

plot(f,P1);
title('Single-Sided Amplitude Spectrum');
xlabel('f (Hz)');
ylabel('|P1(f)|');
axis([0 500 0 5]);
end

