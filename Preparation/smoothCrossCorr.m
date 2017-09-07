function ys = smoothCrossCorr(Xa,Xb,varargin)
%SMOOTHCROSSCORR Calculates the smoothed crosscorrelation of Xa and Xb
%
%   Xa and Xb are quadrature measurements and assumed to be already shaped
%   into piezo-segments (e.g. by prepare3ChData)
%
% Optional Arguments:
%   'Type': You can use a moving average with 'moving' or a cubic smoothing
%   spline with 'spline'. Spline leaves you with more points afterwards,
%   moving average is faster to calculate.

%% Validate and parse input arguments
p = inputParser;
defaultType = 'spline';
defaultParam = 1e-15;
addParameter(p,'Type',defaultType,@isstr);
addParameter(p,'Param',defaultParam,@isvector);
parse(p,varargin{:});
c = struct2cell(p.Results);
[param,type] = c{:};

%% Cross-Correlation of X1 and X2
XProd = Xa.*Xb;

%% Approximate each piezo-segment with the chosen smoothing method
[nPulses,nPieces,nSegments] = size(XProd);

switch type
    % spline is superior, because it doesn't lose points on the boundaries
    case 'spline'
        % Approximate each piezo-segment with a cubic smoothing spline,
        % smoothing parameter given by param should be between 1e-10 and
        % 1e-14. 1e-10 was used previously, 1e-13 yields good results.
        x = 1:(nPulses*nPieces);
        y = reshape(XProd,[nPulses*nPieces nSegments]);
        ys = transpose(csaps(x,y',param,x));
    case 'moving'
        y = reshape(XProd,[nPulses*nPieces nSegments]);
        [ys,~] = moving_average(y,param,1);
        %moving_average by Carlos Vargas.
        %param gives the semilength of the smoothing window, ...
        %i.e. the window has the length 2*param +1. 
        
        % Eliminate boundaries which are averaged over less than 2*param+1
        % points
        ys(1:param,:)=NaN(param,nSegments);
        ys(end-param+1:end,:)=NaN(param,nSegments);
end

end

