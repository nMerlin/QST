function [ I ] = cavityThermalBath()
%Simulates a single mode in a cavity with a surrounding bath

% Parameters
beta = 1;
Omega = 1;
nB = 1; % 1/(exp(beta*Omega)-1) Bose-Einstein bath occupation
Gamma = 0.1; % Coupling constant
nMax = 10; % Maximum matrix dimension
nVector = (-1:nMax+1)';
rho = zeros(nMax+3,1);
rho(nMax+2) = 1;
h = 0.1;
t = 0:h:100;
I = zeros(length(t),1);
rhoSaved = zeros(nMax+3,length(t));
k1 = zeros(nMax+3,1);
k2 = k1;
k3 = k1;
k4 = k1;

% Solve with Runge-Kutta
for j = 1:length(t)
    I(j) = nVector' * rho;
    rhoSaved(:,j)=rho;
    k1(2:end-1) = h * rhoPoint(nVector, rho, Gamma, nB);
    k2(2:end-1) = h * rhoPoint(nVector, rho + 0.5 * k1, Gamma, nB);
    k3(2:end-1) = h * rhoPoint(nVector, rho + 0.5 * k2, Gamma, nB);
    k4(2:end-1) = h * rhoPoint(nVector, rho + k3, Gamma, nB);
    rho = rho + (k1 + 2*k2 + 2*k3 + k4)/6;
    rho = rho / sum(rho);
end

plot(rhoSaved(:,:)');
end

function deriv = rhoPoint(nVector, rho, Gamma, nB)
    deriv = Gamma .* (...
        nB * nVector(2:end-1) .* rho(1:end-2) - ...
        (nVector(2:end-1)+(2*nVector(2:end-1)+1)*nB).*rho(2:end-1) + ...
        (1+nB)*nVector(3:end).*rho(3:end) );
end
