function [g22Ch,times,ada1,ada2] = g22Ch(X1, X2, X12, nResolution1, nResolution2, nResolution12,varargin)
%g22Ch Computes the g2 correlation between two channels channel 1 and
%channel 2. 
%
% Usage:
%    [g22Ch,times,ada1,ada2] = g22Ch(X1, X2, X12, nResolution1, nResolution2, nResolution12,varargin)
%
% Input Arguments:
%   X1, X2 - Quadratures from a continuous quantum state measurement (equal time
%       spacing between all points), of two different channels
%   NRESOLUTION - Number of Quadratures used to create a single data point
%           (time resolution). It is recommended to use  NRESOLUTION =
%           length(X) in order to have enough values where all phases
%           between the two channels are averaged.
%   X12 - Should be the sum of the two channels quadratures: X12 = X1(:)+X2(:);
%   or a uniformly sampled sum with respect to the phase e.g. 
%   X12 = seriesUniformSampling(X1(:)+X2(:),theta12(:),'NBins',100);
%
% Optional Input Arguments:
%   'SampleRate': Sample rate of the recorded quadratures in MHz. Necessary
%       to compute the time-axis.
%
% Output Arguments:
%   g2: g2(0) values
%   ada1, ada2: number of photons in each channel
%   times: recording times in seconds

%% Validate and parse input arguments
p = inputParser;
defaultSampleRate = 75.4; % MHz
addParameter(p,'SampleRate',defaultSampleRate,@isnumeric);
parse(p,varargin{:});
c = struct2cell(p.Results);
[samplerate] = c{:};

%% Reshaping X1 according to NRESOLUTION1
X1 = X1(:);
nSegments = floor(length(X1)/nResolution1);
X1 = X1(1:nSegments*nResolution1);
X1 = reshape(X1,[nResolution1 nSegments]); % Consecutive values in columns

%% Piecewise photon number channel 1 <a^+ a>
X1 = X1 - mean(mean(X1));
ada1 = mean(X1.^2)-0.5;

%% Reshaping X2 according to NRESOLUTION2
X2 = X2(:);
nSegments = floor(length(X2)/nResolution2);
X2 = X2(1:nSegments*nResolution2);
X2 = reshape(X2,[nResolution2 nSegments]); % Consecutive values in columns

%% Piecewise photon number channel 2 <a^+ a>
X2 = X2 - mean(mean(X2));
ada2 = mean(X2.^2)-0.5;

%% Reshaping X12 according to NRESOLUTION12
X12 = X12(:);
nSegments = floor(length(X12)/nResolution12);
X12 = X12(1:nSegments*nResolution12);
X12 = reshape(X12,[nResolution12 nSegments]); % Consecutive values in columns


%% Numerator
Num = 0.25 + (mean(X12.^4) - mean(X1.^4) - mean(X2.^4) - 3*mean(X1.^2) - 3* mean(X2.^2)) /6;

%% Get g2 between X1 and X2
g22Ch = Num./(ada1.*ada2);
g22Ch = g22Ch';

%% Time axis
times = (0.5:1:length(g22Ch))*1/samplerate*nResolution1/1000000; % seconds

end

