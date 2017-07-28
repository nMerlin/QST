function ys = smoothMv(y,span)
%SMOOTHCROSSCORR Calculates the smoothed crosscorrelation of Xa and Xb
%
%   Xa and Xb are quadrature measurements and assumed to be already shaped
%   into piezo-segments (e.g. by prepare3ChData)

%% Handle optional input arguments and default values

%y is of format [nPoints nSegments].
ys = zeros(size(y));

[nPoints, ~] = size(y);
for n= 1:nPoints
    
    if(n < (span + 1)/2)
        ys(n,:) = mean(y(1:2*n-1,:));
    elseif (n > nPoints - (span + 1)/2 )
        ys(n,:)= mean(y(2*n-nPoints:nPoints,:));    
    else
       ys(n,:) = mean( y(n-(span-1)/2:n+(span+1)/2,:) );
    end
        
    
end



end