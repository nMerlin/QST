function [ simData ] = simSampling(pulse,samplerate,reprate,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Validate and parse input arguments
p = inputParser;
defaultNewSamplerate = samplerate;
defaultNPulses = 1000;
defaultPlot = 'hide';
addParameter(p,'NPulses',defaultNPulses,@isnumeric);
addParameter(p,'NewSamplerate',defaultNewSamplerate,@isnumeric);
addParameter(p,'Plot',defaultPlot,@isstr);
parse(p,varargin{:});
c = struct2cell(p.Results);
[newSamplerate,nPulses,plotarg] = c{:};

%% Interpolate given pulse
y = double(pulse);
x = (0:length(y)-1)*1/samplerate;
sx = linspace(0,length(y)-1,100000).*0.2e-9;
sp = csaps(x,y,1);

if strcmp(plotarg,'show')
    plot(x,y,'o'); hold on;
    plot(sx,fnval(sp,sx));
    hold off;
end

%% Simulate sampled pulses
nSamples = floor(x(end)*newSamplerate);
simData = zeros(nSamples,nPulses);
samplevector = (0:nSamples-1)*1/newSamplerate;
for k=1:nPulses
    sampleoffset = mod((k-1)*1/reprate,1/newSamplerate);
    simData(:,k) = fnval(sp,samplevector+sampleoffset);
end
simData = round(simData);

end

