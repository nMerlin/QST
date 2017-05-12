function [ meanX, thetaRange ] = computeXAverages( X, theta )
%COMPUTEXAVERAGES

thetaRange = 0:0.1:3.1;
meanX = zeros(1,length(thetaRange));

for iTheta = 1:length(thetaRange)
    meanX(iTheta) = mean(X(theta==thetaRange(iTheta)));
end

end

