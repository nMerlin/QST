function [solution] = solveCorrProblem(X1,X2,X3)
%UNTITLED2 Summary of this function goes here
%   Comments:
%       This function only works at midnight and full moon on a witch's
%       mountain in a circle which is draw with the blood of a maiden.

nPoints = 900;
[nRows1,nColumns] = size(X1);
[nRows2,~] = size(X2);
[nRows3,~] = size(X3);

for shift = -10:10
    [X1long,X3long] = deal(zeros(nPoints,nColumns));
    X1start = ceil((nRows1-nPoints)/2)+shift;
    X3start = ceil((nRows2-nPoints)/2);
    for iColumn = 1:nColumns
        X1long(:,iColumn) = X1(X1start:X1start+nPoints-1,iColumn);
        X3long(:,iColumn) = X3(X3start:X3start+nPoints-1,iColumn);
    end
    XProd = X1long.*X3long;
    plot(XProd(:));
    w = waitforbuttonpress;
end

end

