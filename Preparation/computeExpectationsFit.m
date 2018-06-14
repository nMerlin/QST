function [delX,varX,meanN,nSegments] = computeExpectationsFit(X,theta,varargin)
%COMPUTEEXPECTATIONSFIT Computes <Q>, <P>, <Q^2>, <P^2> and uncertainties,
%using a sine fit. X and THETA are computed quadratures and phase values
%that are not discretized in THETA.
%
% Output Arguments:
%   delX: 
%
% Optional Input Arguments:
%   'Plot': Set it to 'show' if you want a graphical output. Default is
%       'hide' with no graphical output.
%   'Norm': It's defined by the commutator relation between the two
%       quadrature operators: [q,p] = i*2*Norm
%       Typical values are 1/sqrt(2), 1/2 or 1. Default is 1/sqrt(2).

%% Validate and parse input arguments
p = inputParser;
defaultNorm = 1/sqrt(2);
addParameter(p,'Norm',defaultNorm,@isnumeric);
defaultPlotOpt = 'hide';
addParameter(p,'Plot',defaultPlotOpt,@isstr);
parse(p,varargin{:});
c = struct2cell(p.Results);
[Norm,plotOpt] = c{:};

%% Loop over number of piezo segments
[~,nSegments] = size(X);
[meanN,varX,delX] = deal(zeros(nSegments,1));
for iSegment = 1 : nSegments
    
    % Fit a sine function to model coherent offset
    segmentTheta = theta(:,iSegment);
    segmentX = X(:,iSegment);
    [fitParams,fitFunction] = fitSinusoidal(segmentTheta,segmentX, ...
        'Plot',plotOpt);
    fitValues = fitFunction(fitParams,segmentTheta);
    
    % Compute variance and standard deviation
    fluctuations = segmentX - fitValues;
    varX(iSegment) = var(fluctuations);
    delX(iSegment) = std(fluctuations);

    % Coherent state photon number
    amplitude = fitParams(1);
    meanN(iSegment) = (amplitude^2)/(4*Norm^2);
 
end % iSegment 

end % computeExpectations

