function [ rho ] = computeDensityMatrix( X, theta )
%COMPUTEDENSITYMATRIX Summary of this function goes here
%
%   X and THETA must ...
%       - be 1D-Arrays
%       - have their NaN values at the same places.

% Stripping off NaN values
dispstat('','init');
dispstat('Strip X and THETA off their NaN values.','timestamp','keepthis');
X = X(~isnan(X));
theta = theta(~isnan(theta));
assert(length(X)==length(theta), ['X and THETA should have the same ' ...
    'length after stripping all NaN values!']);

MAX_FOCK_STATE = 30;
N_ITERATIONS = 30;

PI1D = computeProjector1D( X, theta, MAX_FOCK_STATE);
nX = length(X);
prob = zeros(nX, 1);
nextRho = normalize(ones(MAX_FOCK_STATE+1,MAX_FOCK_STATE+1));
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
end

end

