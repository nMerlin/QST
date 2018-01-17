function Xcomp = correlationCompensation(X, varargin)
% The mystical correlation compensator.
%X : Input quadratures. Their offset should be removed already.
%Xcomp : Output quadratures without contamination by precedent pulses. 
% The idea comes from R.Kumar et al: Versatile wideband balanced
% detector..., Optics Communications 285 (2012) 5259-5267

%% Validate and parse input arguments
p = inputParser;
defaultNormalize = 'interior';
addParameter(p,'Normalize',defaultNormalize);
parse(p,varargin{:});
c = struct2cell(p.Results);
[normalize] = c{:};

%% Computation for each channel

[~,~, channels] = size(X);
    
Xcomp = zeros(size(X));
for channel = 1:channels
    Xch = X(:,:,channel)';
    Xcomp(:,:,channel) =  mgs(Xch, normalize)';
end

end

function Q =  mgs(X, normalize)  %The idea comes from modified Gram Schmidt Orthogonalisation
    [n,m] = size(X);
    Q = zeros(n,m);
    N = zeros(n,m);
    for k = 1:m
        Q(:,k) = X(:,k); 
        N(:,k) = X(:,k);
        if strcmp(normalize,'exterior')
        %normalize with Q(k)^2
            R = Q(:,1:k-1)'*Q(:,k) /(Q(:,k)'*Q(:,k));         
            for u = 1:k-1
               Q(:,k) = Q(:,k) - R(u)*Q(:,u);
            end
        end
        
        if strcmp(normalize,'interior')
        %normalize with Q(u)^2
            R = Q(:,1:k-1)'*Q(:,k) ;         
            for u = 1:k-1
               Q(:,k) = Q(:,k) - R(u)*N(:,u);
            end
            N(:,k) = Q(:,k)/(Q(:,k)'*Q(:,k));
        end
        
    end
end
