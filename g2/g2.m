function [g2vec, ada, times] = g2(X, nResolution, varargin)
%G2 Creates a plot showing the g2(0) behavior over time
%
% Usage:
%   [g2,ada,times] = g2(X,nResolution,varargin);
%
% Input Arguments:
%   X - Quadratures from a continuous quantum state measurement (equal time
%       spacing between all points)
%   NRESOLUTION - Number of Quadratures used to create a single data point
%           (time resolution)
%
% Optional Input Arguments:
%   'SampleRate': Sample rate of the recorded quadratures in MHz. Necessary
%       to compute the time-axis.
%
% Output Arguments:
%   g2: g2(0) values
%   ada: number of photons
%   times: recording times in seconds

%% Validate and parse input arguments
p = inputParser;
defaultSampleRate = 75.4; % MHz
addParameter(p,'SampleRate',defaultSampleRate,@isnumeric);
parse(p,varargin{:});
c = struct2cell(p.Results);
[samplerate] = c{:};

%% Reshaping X according to NRESOLUTION
X = X(:);
nSegments = floor(length(X)/nResolution);
X = X(1:nSegments*nResolution);
X = reshape(X,[nResolution nSegments]); % Consecutive values in columns

%% Piecewise photon number <a^+ a>
X = X - mean(mean(X));
ada = mean(X.^2)-0.5;
%ada(ada<0) = 0;

%% Piecewise <a^+ a^+ a a>, in the following adadaa, and g2(t,t)
adadaa = 2/3*mean(X.^4)-2*ada-0.5;
g2vec = adadaa./ada.^2;
g2vec = g2vec';

%% Time axis
times = (0.5:1:length(g2vec))*1/samplerate*nResolution/1000000; % seconds

end

