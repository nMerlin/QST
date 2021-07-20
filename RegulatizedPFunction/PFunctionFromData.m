function [P,sigmaP,QuadVals] = PFunctionFromData(Xdata,thetaData,directory,xGrid,phiGrid)
%calculates the regularized P function from the measured quadratures XData and 
% phases thetaData with help of the pattern functions in directory.
% see Jans notes, equ. 17
%sigmaP is the standard deviation via equ. 18 
%you can get the xGrid and phiGrid from the filenames of the archived
%pattern funcions via 
% [xGrid,phiGrid] = getGridsFromFilenames(directory);

%get x and phi stepsize
xStep = min(diff(xGrid));
phiStep = min(diff(phiGrid));

%get the new coordinate grid from folder name  
maxQuadToken = regexpi(directory,'maxQuad-([-0123456789.]*)-','tokens');
maxQuad = str2double(cell2mat(maxQuadToken{1}));
ResToken = regexpi(directory,'Resolution-([-0123456789.]*)-','tokens');
Resolution = str2double(cell2mat(ResToken{1}));
QuadVals=-abs(maxQuad):Resolution:abs(maxQuad);

% compute the P function and its uncertainty 
summe = 0;
sum2 = 0;
for phi = phiGrid
    for x = xGrid 
        load(strcat(directory,'\x-',num2str(x),'-phi-',num2str(phi),'.mat'),'pattern');
        binX = [x-xStep/2, x+xStep/2];
        binPhi = [phi-phiStep/2, phi+phiStep/2];
        N = length(Xdata(Xdata > min(binX) & Xdata < max(binX) & thetaData > min(binPhi) & thetaData < max(binPhi)));
        N2 = length(Xdata(thetaData > min(binPhi) & thetaData < max(binPhi)));
        weight = N/N2;
        summe = summe + weight * pattern;
        sum2 = sum2 + weight * pattern.^2;
    end
end
P = summe * pi / length(phiGrid);
P2 = sum2 * pi^2 / length(phiGrid);
sumP = sum(sum(P));
P=P./sumP;

Ntotal = length(Xdata);
sigmaP = sqrt((P2 - P.^2)/(Ntotal-1));
sigmaP=sigmaP./sumP;
end


