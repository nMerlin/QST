function out = exponentialCDF(in,nAv,varargin)
%EXPONENTIALCDF Cumulative distribution function (and inverse) of
%exponential distribution for average number of photons nAv

%% Validate and parse input arguments
p = inputParser;
defaultInverse = false;
addParameter(p,'Inverse',defaultInverse,@islogical);
parse(p,varargin{:});
c = struct2cell(p.Results);
[inverse] = c{:};

%% Compute bose-einstein cumulative distribution function or inverse
lambda = 1./nAv;
if inverse
    out = -1/lambda.*log(1-in);
else
    out = 1 - exp(-lambda.*in);
end

end