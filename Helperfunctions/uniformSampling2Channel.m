function [X1Start,X2Start,theta1Start,theta2Start] = uniformSampling2Channel(X1,X2,theta1,theta2,iterN )
%UNTITLED2 Summary of this function goes here
%  This function ensures both X1 and X2 have all phases uniformed AND their relative phase uniformed.
% X1, X2: quadratures from homodyne measurement.
% theta1, theta2: there phases that were computed with computePhase.m
% iterN: number of iterations of the uniforming process. 

X1Start = X1(:);
X2Start = X2(:);
theta1Start = theta1(:);
theta2Start = theta2(:);
theta12Start = mod(theta2(:)-theta1(:),2*pi);

for i = 1: iterN
    %uniforming of relative phase
    [~,~,uniformIndicesMidX12] = seriesUniformSamplingIndex(X1Start+X2Start,theta12Start,'NBins',300);
    X1Start = X1Start(uniformIndicesMidX12);
    X2Start = X2Start(uniformIndicesMidX12);
    theta1Start = theta1Start(uniformIndicesMidX12);
    theta2Start = theta2Start(uniformIndicesMidX12);

    %another uniforming of channel 1 phase 
    [X1Start,theta1Start,indicesX1End] = seriesUniformSamplingIndex(X1Start,theta1Start,'NBins',300);
    X2Start = X2Start(indicesX1End);
    theta2Start = theta2Start(indicesX1End);

    %another uniforming of channel 2 phase 
    [X2Start,theta2Start, indicesX2End] = seriesUniformSamplingIndex(X2Start,theta2Start,'NBins',300);
    X1Start = X1Start(indicesX2End);
    theta1Start = theta1Start(indicesX2End);

    theta12Start = mod(theta2Start-theta1Start,2*pi);
end

end

