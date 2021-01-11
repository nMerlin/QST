function [ X, theta ] = discretizeTheta( X, theta, nIntervals )
%DISCRETIZETHETA Divides the range of theta values into discrete intervals

[nRows, nSegments] = size(theta);
XOut = NaN(ceil(3*nRows/nIntervals), nIntervals, nSegments);
thetaOut = XOut;
maxCount = 0; % Keeps track of longest interval
for iSeg = 1 : nSegments
    % Remove NaN entries    
    thetaSegment = theta(~(isnan(theta(:,iSeg))) & ~(isnan(X(:,iSeg))),iSeg);
    XSegment = X(~(isnan(theta(:,iSeg))) & ~(isnan(X(:,iSeg))),iSeg);
    
    % Sorting theta and X into intervals with MatLab functions 
    [N,~,bin] = histcounts(thetaSegment,nIntervals);
    [~,I] = sort(bin);
    thetaSegment = thetaSegment(I);
    XSegment = XSegment(I);
    
    for iInterval = 1 : nIntervals
        start = 1+sum(N(1:iInterval-1));
        stop = start+N(iInterval)-1;
        thetaOut(1:N(iInterval),iInterval,iSeg) = thetaSegment(start:stop);
        XOut(1:N(iInterval),iInterval,iSeg) = XSegment(start:stop);
    end
    maxCount = max(max(N), maxCount);
end

% Clean output variables from unnecessary NaN values
X = XOut(1:maxCount,:,:);
X(X==0) = NaN;
theta = thetaOut(1:maxCount,:,:);
theta(theta==0) = NaN;

end

