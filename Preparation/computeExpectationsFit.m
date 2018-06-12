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

Norm = 1/sqrt(2);
%Norm ist the factor in the relation between the quadratures and the ladder
%operators: q = Norm*(a^{+} + a), p = Norm*i*(a^{+} - a)
%typical values are 1/sqrt(2) or 1/2. 

%% Validate and parse input arguments
p = inputParser;
defaultPlotOpt = 'hide';
addParameter(p,'Plot',defaultPlotOpt,@isstr);
parse(p,varargin{:});
c = struct2cell(p.Results);
[plotOpt] = c{:};

%% Loop over number of piezo segments
[~,nSegments] = size(X);
[meanN,varX,delX] = deal(zeros(nSegments,1));
for iSegment = 1 : nSegments
    % Compute the Fit
    segmentTheta = theta(:,iSegment);
    segmentX = X(:,iSegment);
    
    [fitParams,fitFunction,~,~] = fitSinus(segmentTheta,segmentX);
    functionValues=fitFunction(fitParams,segmentTheta);
    X0 = fitParams(1);
    
    deviation = abs(segmentX - functionValues);
    
    %compute the variance and std interval-wise
    
    %compute the total variance and std
    varX(iSegment) = sum(deviation.^2)/(length(deviation)-1) ;
    delX(iSegment) = sqrt(varX(iSegment));

    % Coherent state photon number
    meanN(iSegment) = (X0^2)/2;
 
end % iSegment 

end % computeExpectations

