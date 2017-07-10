function ys = smoothCrossCorr(Xa,Xb,varargin)
%SMOOTHCROSSCORR Calculates the smoothed crosscorrelation of Xa and Xb
%
%   Xa and Xb are quadrature measurements and assumed to be already shaped
%   into piezo-segments (e.g. by prepare3ChData)

%% Handle optional input arguments and default values
nVarargin = length(varargin);
optArgs = {0.0000000001};
optArgs(1:nVarargin) = varargin;
[P_SMOOTH] = optArgs{:}; % Smoothing parameter

%% Cross-Correlation of X1 and X2
XProd = Xa.*Xb;

%% Approximate each piezo-segment with a cubic smoothing spline
[nPulses,nPieces,nSegments] = size(XProd);
x = 1:(nPulses*nPieces);
y = reshape(XProd,[nPulses*nPieces nSegments]);
ys = transpose(csaps(x,y',P_SMOOTH,x));
ys = ys(:);

end

