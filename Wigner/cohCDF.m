function out = cohCDF(in,theta,nAv,varargin)
%COHCDF Cumulative distribution function of integral projections from
%arbitrary coherent Wigner function.
%
%   Input Arguments:
%       in: Input value or vector to apply the function on.
%       theta: Phase-space angle to the projection axis in degrees.
%       nAv: Average number of photons in the coherent state.
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

%% Compute coherent cumulative distribution function or inverse
if inverse
    out = sqrt(2*nAv).*cosd(theta)+erfinv(-1+2*in);
else
    out = 1/2*(1 + erf(in-sqrt(2*nAv).*cosd(theta)));
end

end

