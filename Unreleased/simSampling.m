function [ simData ] = simSampling(pulse,reprate,oldsamplerate,newsamplerate,varargin)
%SIMSAMPLING Simulates the sampling of a real pulse train.
%   
%   Output Arguments:
%     SIMDATA: Array containing consecutive sampled pulses
%
%   Input Arguments:
%     PULSE: Example pulse sampled with OLDSAMPLERATE. This data will be
%       interpolated.
%     REPRATE: Repetition rate of the pulsed laser. Necessary to compute
%       shifts due to sampling.
%     NEWSAMPLERATE: This samplerate will be used for SIMDATA.
%
%   Optional Input Arguments:
%     'NPulses': Number of pulses to simulate.
%     'Plot': With 'show' the function plots the interpolated PULSE.

%% Validate and parse input arguments
p = inputParser;
defaultNPulses = 1000;
defaultPlot = 'hide';
addParameter(p,'NPulses',defaultNPulses,@isnumeric);
addParameter(p,'Plot',defaultPlot,@isstr);
parse(p,varargin{:});
c = struct2cell(p.Results);
[nPulses,plotarg] = c{:};

%% Interpolate given pulse
y = double(pulse);
x = (0:length(y)-1)*1/oldsamplerate;
sx = linspace(0,length(y)-1,100000).*0.2e-9;
sp = csaps(x,y,1);

if strcmp(plotarg,'show')
    plot(x,y,'o'); hold on;
    plot(sx,fnval(sp,sx));
    hold off;
end

%% Simulate sampled pulses
nSamples = floor(x(end)*newsamplerate);
simData = zeros(nSamples,nPulses);
samplevector = (0:nSamples-1)*1/newsamplerate;
for k=1:nPulses
    sampleoffset = mod((k-1)*1/reprate,1/newsamplerate);
    simData(:,k) = fnval(sp,samplevector+sampleoffset);
end
simData = round(simData);

end

