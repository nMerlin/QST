function [ rho ] = computeDensityMatrix( X, theta )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

MAX_FOCK_STATE = 20;
N_ITERATIONS = 20;

PI1D = computeProjector1D( X, theta, MAX_FOCK_STATE);
nX = length(X);
prob = zeros(nX, 1);
nextRho = normalize(ones(MAX_FOCK_STATE+1,MAX_FOCK_STATE+1));
for iRho = 1:N_ITERATIONS
    % Iteration operator R = sum_i(Projector(theta_i,x_i)/prob_thetai(x_i))
    % Iteration: rho_k+1 = norm(R rho_k R)
    iRho
    rho = nextRho;
    R = zeros(MAX_FOCK_STATE+1,MAX_FOCK_STATE+1);
    for i = 1:nX
        prob(i) = PI1D(:,i)' * rho * PI1D(:,i);
        R = R + (PI1D(:,i)*PI1D(:,i)')/(nX *prob(i));
    end
    nextRho = R * rho * R; % ITERATION step
    nextRho = normalize(nextRho); % normalization
end

end

