function out = cohCDF(in,theta,nAv,varargin)
%COHCDF Cumulative distribution function of integral projections from
%arbitrary coherent Wigner functions

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

