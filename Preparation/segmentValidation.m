function [X1,X2,X3,theta] = segmentValidation(X1,X2,X3,theta,varargin)
%SEGMENTVALIDATION manual selection of correctly reconstructed phases
%
% Input Parameters:
%   X1,X2,X3 - quadrature data in the format [nPoints,nRecords,nSegments]
%   theta - reconstructed phase with cross correlation X1.*X3
%
% Output Paramters:
%   X1,X2,X3,theta - only the manually accepted parts of the corresponding
%       input arrays are returned

%% Validate and parse input arguments
p = inputParser;
defaultFilename = '';
addParameter(p,'Filename',defaultFilename,@isstr);
parse(p,varargin{:});

end

