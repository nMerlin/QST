function [rho,history] = computeDensityMatrix( X, theta, varargin )
%COMPUTEDENSITYMATRIX Summary of this function goes here
%
%   X and THETA must ...
%       - be 1D-Arrays
%       - have their NaN values at the same places.

%% Constants
MAX_FOCK_STATE = 20;
N_ITERATIONS = 100;

%% Validate and parse input arguments
p = inputParser;
defaultIterations = N_ITERATIONS;
defaultRho = normalize(ones(MAX_FOCK_STATE+1,MAX_FOCK_STATE+1));
addParameter(p,'Rho',defaultRho,@ismatrix);
addParameter(p,'Iterations',defaultIterations,@isnumeric);
parse(p,varargin{:});
c = struct2cell(p.Results);
[nextRho,N_ITERATIONS] = c{:};

% Stripping off NaN values
dispstat('','init');
dispstat('Strip X and THETA off their NaN values.','timestamp','keepthis');
X = X(~isnan(X));
theta = theta(~isnan(theta));
assert(length(X)==length(theta), ['X and THETA should have the same ' ...
    'length after stripping all NaN values!']);

%% History
% The _history_ variable will contain all intermediately calculated density
% matrices. It is only saved when the output is required.
if nargout>1
    history = zeros(MAX_FOCK_STATE+1,MAX_FOCK_STATE+1,N_ITERATIONS);
end

%% Iteration
PI1D = computeProjector1D( X, theta, MAX_FOCK_STATE);
nX = length(X);
prob = zeros(nX, 1);
for iRho = 1:N_ITERATIONS
    % Iteration operator R = sum_i(Projector(theta_i,x_i)/prob_thetai(x_i))
    % Iteration: rho_k+1 = norm(R rho_k R)
    dispstat(['Iteration: ' num2str(iRho) ' of ' ...
        num2str(N_ITERATIONS) '.'], 'timestamp');
    rho = nextRho;
    R = zeros(MAX_FOCK_STATE+1,MAX_FOCK_STATE+1);
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
end

end

