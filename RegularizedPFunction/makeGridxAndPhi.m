function [xGrid,phiGrid] = makeGridxAndPhi(maxX,xStep,phiStep)
    xGrid = -abs(maxX):xStep:abs(maxX);
    phiGrid = 0:phiStep:2*pi;
end