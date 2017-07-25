function [h] = plotRho(rho,varargin)
%PLOTRHO Plots the real part of the complex valued density matrix RHO
% 
% Output Arguments:
%   H - handle of the generated figure

%% Validate and parse input arguments
p = inputParser;
defaultHandle = gca;
addParameter(p,'Handle',defaultHandle,@ishandle);
parse(p,varargin{:});
c = struct2cell(p.Results);
[axes_handle] = c{:};

h = bar3(axes_handle,real(rho));

end

