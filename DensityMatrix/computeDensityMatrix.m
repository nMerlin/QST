function [rho,history] = computeDensityMatrix( X, theta, varargin )
%COMPUTEDENSITYMATRIX Summary of this function goes here
%
%   X and THETA must ...
%       - be 1D-Arrays
%       - have their NaN values at the same places.
%
% Optional Input Arguments:
%   'Iterations': Default is 100. How many iterations take place in a loop.
%   'Threshold': EXPERIMENTAL. Default is 0. Stop loop,
%       when |rho - nextrho| < threshold.

%% Validate and parse input arguments
p = inputParser;
defaultIterations = 100;
addParameter(p,'Iterations',defaultIterations,@isnumeric);
defaultMaxFockState = 30;
addParameter(p,'MaxFockState',defaultMaxFockState,@isnumeric);
defaultRho = [];
addParameter(p,'Rho',defaultRho,@ismatrix);
defaultThreshold = 0;
addParameter(p,'Threshold',defaultThreshold,@isnumeric);
parse(p,varargin{:});
c = struct2cell(p.Results);
[N_ITERATIONS,maxFockState,nextRho,threshold] = c{:};

if isempty(nextRho)
    nextRho = normalize(ones(maxFockState+1,maxFockState+1));
end

% Stripping off NaN values
dispstat('','init');
dispstat('Strip X and THETA off their NaN values.','timestamp','keepthis');
if sum(isnan(X))~=sum(isnan(theta))
    warning(['X and THETA should have the same ',...
        'length after stripping all NaN values!']);
    X = X(~isnan(theta));
else
    X = X(~isnan(X));
end
theta = theta(~isnan(theta));

%% History
% The _history_ variable will contain all calculated density
% matrices. It is only saved when the output is required.
if nargout>1
    history = zeros(maxFockState+1,maxFockState+1,N_ITERATIONS);
end

%% Iteration
PI1D = computeProjector1D( X, theta, maxFockState);
nX = length(X);
prob = zeros(nX, 1);
for iRho = 1:N_ITERATIONS
    % Iteration operator R = sum_i(Projector(theta_i,x_i)/prob_thetai(x_i))
    % Iteration: rho_k+1 = norm(R rho_k R)
    dispstat(['Iteration: ' num2str(iRho) ' of ' ...
        num2str(N_ITERATIONS) '.'], 'timestamp');
    rho = nextRho;
    R = zeros(maxFockState+1,maxFockState+1);
    for i = 1:nX
        prob(i) = PI1D(:,i)' * rho * PI1D(:,i);
        R = R + (PI1D(:,i)*PI1D(:,i)')/(nX*prob(i));
    end
    nextRho = R * rho * R; % ITERATION step
    nextRho = normalize(nextRho); % normalization
    
    % Update history output
    if nargout>1
        history(:,:,iRho) = nextRho;
    end
    
    % Stop when change is smaller than threshold
    if sum(sum(abs(rho-nextRho))) < threshold
        rho = nextRho;
        dispstat(['Iteration: ',num2str(iRho),' of ', ...
            num2str(N_ITERATIONS),'. Threshold reached!'], ...
            'timestamp','keepthis');
        break;
    end
end

end

