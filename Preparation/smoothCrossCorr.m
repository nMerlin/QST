function ys = smoothCrossCorr(Xa,Xb,varargin)
%SMOOTHCROSSCORR Calculates the smoothed crosscorrelation of Xa and Xb
%
%   Xa and Xb are quadrature measurements and assumed to be already shaped
%   into piezo-segments (e.g. by prepare3ChData)
%   The smoothing method can be chosen by the name-value-pair 'Type','chosen Type'.

%% Validate and parse input arguments
p = inputParser;
defaultType = 'moving';
defaultParam = 5000;
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
    case 'moving'
        y = reshape(XProd,[nPulses*nPieces*nSegments 1]);
        [ys,~] = moving_average(y,param,1); %smoothing with moving average...
        %by Carlos Vargas.
        %param gives the semilength of the smoothing window, ...
        %i.e. the window has the length 2*param +1. 
        
        ys(1:0.1*param)=moving_average(y(1:0.1*param),0.01*param,1);
        
        ys = reshape(ys,[nPulses*nPieces nSegments]);
    case 'spline'
        x = 1:(nPulses*nPieces);
        y = reshape(XProd,[nPulses*nPieces nSegments]);
        ys = transpose(csaps(x,y',param,x)); %Approximate each piezo-segment ..
        %with a cubic smoothing spline, smoothing parameter given by param
        %should be between 1e-10 and 1e-14. 1e-10 was used previously, 1e-13
        %yields good results.
end

end

