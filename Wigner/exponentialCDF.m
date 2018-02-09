function out = exponentialCDF(in,exp,varargin)
%EXPONENTIALCDF Cumulative distribution function (and inverse) of
%exponential distribution.
%
%   Input Arguments:
%       in: Input value or vector to apply the function on.
%       exp: Expectation value of the exponential distribution.
%
%   Optional Input Arguments:
%       'Inverse': Boolean value indicating whether to use the inverse CDF
%         or the CDF itself. Default is false (not inverting CDF).
%
%   Output Arguments:
%       out: Output value or vector after applying CDF or inverse CDF.

%% Validate and parse input arguments
p = inputParser;
defaultInverse = false;
addParameter(p,'Inverse',defaultInverse,@islogical);
parse(p,varargin{:});
c = struct2cell(p.Results);
[inverse] = c{:};

%% Compute exponential cumulative distribution function or inverse
lambda = 1./exp;
if inverse
    out = -1/lambda.*log(1-in);
else
    out = 1 - exp(-lambda.*in);
end

end