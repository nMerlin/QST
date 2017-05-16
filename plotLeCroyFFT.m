function [] =  plotLeCroyFFT(filename)

cd('raw-data');
M = dlmread(filename, ',',5,0);
f = M(:,1);
fft = M(:,2);
cd('..');

Npeaks = 3;

% Plotting
[pks,locs] = findpeaks(fft,f,'MinPeakProminence',max(abs(fft))/4, 'Npeaks',Npeaks);  %,'MinPeakProminence',max(fft)/10

plot(f,fft);
hold on;
plot(locs,pks+2,'v','markerFaceColor','b','markerEdgeColor','b');
for i = 1:min(Npeaks,length(locs))
    text(locs(i)+10,pks(i),num2str(locs(i)));
end
title('Single-Sided Amplitude Spectrum from LeCroy');
xlabel('f (Hz)');
ylabel('Spectral Power (dB)');
axis([0 500 min(fft) 0]);
grid();

% Saving 
outputFilename = ['LeCroy' filename '.jpg'];
outputFiletype = '-djpeg';
print(outputFilename,outputFiletype);

end
