function [ locs ] = plotFFT( filename )
%PLOTFFT taken from https://de.mathworks.com/help/matlab/ref/fft.html
%   SIGNAL: Recorded Spectraum
%   SAMPLERATE: Sampling Frequency

% Loading
[data8bit,config,~]=load8BitBinary(filename,'dontsave');
originalSamplerate = config.SpectrumCard.Clock.SamplingRate0x28MHz0x29_DBL;
signal = mean(data8bit(:,:,1));
samplerate = round(originalSamplerate*1e6/size(data8bit,1));

% Calculating
T = 1/samplerate;
L = length(signal);
Y = fft(signal);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = samplerate*(0:(L/2))/L;

% Plotting
[pks,locs] = findpeaks(P1,f,'MinPeakProminence',max(P1)/10);
findpeaks(P1,f,'MinPeakProminence',max(P1)/10);
text(locs+10,pks,num2str(locs'));
title('Single-Sided Amplitude Spectrum');
xlabel('f (Hz)');
ylabel('|P1(f)|');
axis([0 500 0 max(P1)]);

% Saving 
outputFilename = ['FFT-' filename '.jpg'];
outputFiletype = '-djpeg';
print(outputFilename,outputFiletype);
end