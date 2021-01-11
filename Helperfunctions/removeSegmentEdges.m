function [X1rem,X2rem,X3rem,thetaRem,n1rem,n2rem,n3rem] = removeSegmentEdges(X1,X2,X3,theta,n1,n2,n3)
%%removes the edges of the piezo segments, where the photon number vectors make jumps
edge = 0.02;
X1rem = rem(X1,edge);
X2rem = rem(X2,edge);
X3rem = rem(X3,edge);
n1rem = rem(n1,edge);
n2rem = rem(n2,edge);
n3rem = rem(n3,edge);
[d,~] = size(theta);
thetaRem = theta(d*edge:d*(1-edge),:);

    function [xrem] = rem(x,edge)
        [a,b,c] = size(x);
        xrem = reshape(x,[a*b,c]);
        xrem = xrem(a*b*edge:a*b*(1-edge),:);
        [a,~]=size(xrem);
        xrem = reshape(xrem,[a,1,c]);
        
    end

end