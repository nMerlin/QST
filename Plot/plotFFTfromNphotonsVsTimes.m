signal = ada;
T = mean(diff(times));
samplerate = 1/T;

% Calculating
L = length(signal);
Y = fft(signal);
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
axis([0 7000 0 0.5]);

% Saving 
% outputFilename = ['FFT-' filename '.jpg'];
% outputFiletype = '-djpeg';