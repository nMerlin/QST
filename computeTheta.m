function [ X, theta ] = computeTheta( X )
%COMPUTETHETA Computes X and THETA ready for the reconstruction algorithm
%
%   X should have the size [nPulses, nRecords, nSegments]

XOld = X;
[nPulses, nRecords, nSegments] = size(XOld);
X = zeros(nPulses * nRecords, nSegments);
theta = X;
for iSeg = 1:nSegments
    % Average over recordings and create fit input x/y data
    yFit = mean(XOld(:,:,iSeg));
    xFit = (1:nRecords) * nPulses - round(nPulses/2);
    xFit(isnan(yFit)) = NaN;
    
    % Somehow the fit only works for the correct x magnitude
    xFitMagnitude = ceil(log10(max(xFit)));
    xFit = xFit / 10^(xFitMagnitude); % scale x-axis for fitting routine
    
    % Fitting
    [fitParams, ~] = fitSinusoidal(xFit, yFit);
    
    X(:,iSeg) = reshape(XOld(:,:,iSeg), nPulses * nRecords,1);
    xTheta = 1 : nPulses * nRecords;
    theta(:,iSeg) = ...
        mod(2 * pi / fitParams(2) * xTheta / 10^(xFitMagnitude) + ...
        2 * pi / fitParams(3) + pi / 2, 2 * pi);
    theta(isnan(X)) = NaN;
end

end

