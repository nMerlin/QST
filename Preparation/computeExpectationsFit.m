function expectations = computeExpectationsFit(X,theta,varargin)
%COMPUTEEXPECTATIONSFIT Computes <Q>, <P>, <Q^2>, <P^2> and more
%expectation values from a quadarature measurement. X and THETA are
%computed quadratures and phase values that are not discretized in THETA.
%
% Output Arguments:
%   expectations: Structure containing the variance 'varX' of X at a fixed
%       phase, the coherent amplitude 'cohAmpl' and the number of photons
%       'cohN' estimated from the coherent amplitude, together with other
%       expectation values requested with the optional input argument
%       'AddOutputs'.
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
[expectations.cohN,expectations.varX, ...
    expectations.cohAmpl] = deal(zeros(nSegments,1));
for iSegment = 1 : nSegments
    % Fit a sine function to model coherent offset
    segmentTheta = theta(:,iSegment);
    segmentX = X(:,iSegment);
    [fitParams,fitFunction] = fitSinusoidal(segmentTheta,segmentX, ...
        'Plot',plotOpt);
    fitValues = fitFunction(fitParams,segmentTheta);
    
    % Compute variance and standard deviation
    fluctuations = segmentX - fitValues;
    expectations.varX(iSegment) = var(fluctuations);

    % Coherent state photon number
    expectations.cohAmpl(iSegment) = fitParams(1);
    expectations.cohN(iSegment) = ...
        (expectations.cohAmpl(iSegment)^2)/(4*Norm^2);
 
end % iSegment

end % computeExpectations

