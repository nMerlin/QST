function [ X, theta ] = discretizeTheta( X, theta, nIntervals )
%DISCRETIZETHETA Divides the range of theta values into discrete intervals

EPSILON = 1e-12;

[nRows, nSegments] = size(theta);
XOut = NaN(round(nRows / nIntervals * 2), nIntervals, nSegments);
thetaOut = XOut;
maxCount = 0; % Keeps track of longest interval
for iSeg = 1 : nSegments
    thetaSegment = theta(:,iSeg);
    
    % Remove NaN entries
    thetaSegment = thetaSegment(~isnan(thetaSegment));
    
    % Compute interval length
    offset = min(thetaSegment) - EPSILON;
    range = max(thetaSegment) - min(thetaSegment) + nIntervals * EPSILON;
    interval = range / nIntervals;
    
    % Sort theta and X into these intervals
    intervalIndices = ceil((thetaSegment - offset) / interval);
    I = ones(nIntervals,1);
    for iPoint = 1 : length(intervalIndices)
        iInterval = intervalIndices(iPoint);
        thetaOut(I(iInterval), iInterval, iSeg) = theta(iPoint, iSeg);
        XOut(I(iInterval), iInterval, iSeg) = X(iPoint, iSeg);
        I(iInterval) = I(iInterval) + 1;
    end
    maxCount = max(max(I), maxCount);
end
maxCount = maxCount - 1; % Correct for +1 after adding any element

% Clean output variables from unnecessary NaN values
X = XOut(1:maxCount,:,:);
theta = thetaOut(1:maxCount,:,:);

end

