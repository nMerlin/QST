function [  ] = piecewiseAnalysis( X )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nPieces = size(X,2);
[nAv,offset] = deal(zeros(nPieces,1));
for iPiece = 1:nPieces
    % Remove offset
    P = X(:,iPiece);
    P = P - mean(mean(X));
    
    % Mean photon number
    nAv(iPiece) = mean(P.^2)-0.5;
    offset(iPiece) = mean(P);
end

xAxis = 1:length(nAv);
plot(xAxis,nAv,xAxis,offset*100);

end

