function [ ] = plotOscilloscope(varargin)
%plots oscilloscope traces and makes fourier trafo

%% Discover *.csv files
filestruct = dir('raw-data/*.csv');
files = {filestruct.name};


%% Iterate through data files
filerange = 1:length(files);
for i = 1:length(filerange)
     C = strsplit(files{i},'.');
     filename = C{1};
     m = csvread(['raw-data/' filename '.csv'],2); 
     tInc = csvread(['raw-data/' filename '.csv'],1,2,[1,2,1,2]); %time increment
     P = m(:,1);  %measured signal
     t = (1:length(P))*tInc;
     
     %plot time signal
     plot(t*1000,P);
     xlabel('time (ms)');
     ylabel('signal voltage');
     graphicsSettings;
     title(filename);
     savefig([filename '-time.fig']);
     print([filename '-time.png'], '-dpng');
     
     %Fouriertrafo 
     plotFFTfromNphotonsVsTimes(P,t,filename,'Axis',[0 5000 0 0.05]);
end
    
end
