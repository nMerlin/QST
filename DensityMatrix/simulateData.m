function [ X, theta ] = simulateData( )
%SIMULATEDATA simulates a quantum state tomography experiment

N_X_PER_ANGLE = 10000;
QUADRATURE_RANGE = -20:0.01:20;
PHASE_RANGE = 0:0.1:pi;

% Density matrix to simulate
rho = zeros(15,15);
rho(1,1) = 0;
rho(2,2) = 0;
rho(3,3) = 0;
rho(5,5) = 0;
rho(10,10) = 1;
rho = rho/trace(rho);

% Discretize X, theta
theta = PHASE_RANGE;
X = zeros(N_X_PER_ANGLE,length(theta));
quadratures = zeros(length(QUADRATURE_RANGE),length(theta));
prob = quadratures;
for iTheta = 1:length(theta)
    % compute probability density for the current theta-value
    quadratures(:,iTheta) = QUADRATURE_RANGE;
    PI1D = computeProjector1D(quadratures(:,iTheta), theta(iTheta), 14);
    for iQuad=1:length(quadratures)
        prob(iQuad,iTheta) = PI1D(:,iQuad)' * rho * PI1D(:,iQuad);
    end
    prob(:,iTheta) = prob(:,iTheta)/sum(prob(:,iTheta));
    
    % compute quadrature values according to this probability distribution
    X(:,iTheta) = distributeRandomly(N_X_PER_ANGLE, ...
        quadratures(:,iTheta), prob(:,iTheta));
end

% build X and theta as column vectors
X = X(:);
theta = ones(N_X_PER_ANGLE,1) * theta;
theta = theta(:);

end

