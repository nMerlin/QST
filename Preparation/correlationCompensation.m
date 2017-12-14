function Xcomp = correlationCompensation(X)
% The mystical correlation compensator.
%X : Input quadratures. Their offset should be removed already.
%Xcomp : Output quadratures without contamination by precedent pulses. 
% The idea comes from R.Kumar et al: Versatile wideband balanced
% detector..., Optics Communications 285 (2012) 5259-5267

%% Computation for each channel

[~, ~, channels] = size(X);
    
Xcomp = zeros(size(X));
for channel = 1:channels
    Xch = X(:,:,channel)';
    Xcomp(:,:,channel) =  mgs(Xch)';        
end

end

function Q =  mgs(X)  %The idea comes from modified Gram Schmidt Orthogonalisation
    [n,m] = size(X);
    Q = zeros(n,m);
    for k = 1:m
        Q(:,k) = X(:,k);   
        R = Q(:,1:k-1)'*Q(:,k) /(Q(:,k)'*Q(:,k));         
        for u = 1:k-1
           Q(:,k) = Q(:,k) - R(u)*Q(:,u);
        end
    end
end
