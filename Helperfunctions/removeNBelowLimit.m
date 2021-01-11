function [X1rem,X2rem,X3rem,thetaRem,iSelect] = removeNBelowLimit(X1,X2,X3,theta,n1,n2,n3,range)
% removes those Xi where the photon number is below a limit
  
[lowerLim1,upperLim1] = limits(n1,range);
[lowerLim2,upperLim2] = limits(n2,range);
[lowerLim3,upperLim3] = limits(n3,range);
    
iSelect = find(n1 >= lowerLim1 & n1 <= upperLim1  ...
     &  n2 >= lowerLim2 & n2 <= upperLim2 ...
    &   n3 >= lowerLim3 & n3 <= upperLim3 );

X1rem = X1(iSelect);
X2rem = X2(iSelect);
X3rem = X3(iSelect);
thetaRem = theta(iSelect);

    function [lowerLim,upperLim] = limits(n,range)
        lowerLim = max(n(:))*(1-range);
        upperLim = max(n(:));
    end


end