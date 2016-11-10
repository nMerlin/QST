function [ PI1D ] = computeProjector1D( X, theta, maxFockState )
%COMPUTEPROJECTOR1D <n|X,theta> for n = 0:maxFockState

PI1D =zeros(maxFockState+1,length(X));
for n=0:maxFockState
    PI1D(n+1,:) =fockstate(n,X) .*exp(1i*n*theta);
end

end

