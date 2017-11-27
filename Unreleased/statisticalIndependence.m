function [results] = statisticalIndependence(filename)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[data8bit,config,~]=load8BitBinary(filename,'dontsave');
X = computeQuadratures(data8bit,config,1);
X1 = X(:,:,1);
X2 = X(:,:,2);
X3 = X(:,:,3);

results = zeros(3);

for delay = 1:3
    results(1,delay) = delayCorr(X1,delay);
    results(2,delay) = delayCorr(X2,delay);
    results(3,delay) = delayCorr(X3,delay);
end

end

function corr = delayCorr(X, delay)
    A = X(1:end-delay);
    B = X(1+delay:end);
    corr = corrcoef(A,B);
    corr = corr(1,2);
end
