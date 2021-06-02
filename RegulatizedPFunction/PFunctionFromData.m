function [P,sigmaP] = PFunctionFromData(Xdata,thetaData,directory,xGrid,phiGrid)
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

sum = 0;
sum2 = 0;
for phi = phiGrid
    for x = xGrid 
        load(strcat(directory,'\x-',num2str(x),'-phi-',num2str(phi),'-Pattern.mat'),'pattern');
        binX = [x-xStep/2, x+xStep/2];
        binPhi = [phi-phiStep/2, phi+phiStep/2];
        u = Xdata > min(binX) & Xdata < max(binX) & thetaData > min(binPhi) & thetaData < max(binPhi);
        N =  length(u(u>0));
        u2 = thetaData > min(binPhi) & thetaData < max(binPhi);
        N2 = length(u2(u2>0));
        weight = N/N2;
        sum = sum + weight * pattern;
        sum2 = sum2 + weight * pattern.^2;
    end
end
P = sum * pi / length(phiGrid);
P2 = sum2 * pi^2 / length(phiGrid);

Ntotal = length(Xdata);
sigmaP = sqrt((P2 - P.^2)/(Ntotal-1));

end


