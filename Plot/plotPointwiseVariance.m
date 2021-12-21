function [  ] = plotPointwiseVariance( data, varargin )
% PLOTPOINTWISEVARIANCE calculates the pointwise variance of DATA, plots
% the pointwise variance data and highlights the maxima that are above the
% mean value. The columns of DATA are the datasets.
%
%   PLOTPOINTWISEVARIANCE(DATA) Plots the pointwise variance of DATA.
%
%   PLOTPOINTWISEVARIANCE(DATA, FILENAME) Plots the pointwise variance of
%   DATA and outputs the resulting graph to FILENAME. Does the filename
%   include an extenstion, this will be used as the filetype. If not, 'png'
%   is the standard filetype.
%   
%   See also: POINTWISEVARIANCE

%% Validate and parse input arguments
p = inputParser;
defaultMinPeakDistance = 10;
addParameter(p,'MinPeakDistance',defaultMinPeakDistance,@isnumeric);
parse(p,varargin{:});
c = struct2cell(p.Results);
[minPeakDist] = c{:};

%% Plot pointwise variance
pointwiseVariance(data,'Plot','show', 'MinPeakDistance',minPeakDist);



end