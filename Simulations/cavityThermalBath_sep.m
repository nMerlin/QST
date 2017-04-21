function [ I ] = cavityThermalBath_sep()
%Simulates a single mode in a cavity with a surrounding bath

% Parameters
outputFilename = 'CavityThermalBath.jpg';
outputFiletype = '-djpeg';
outputFilename2 = 'Intensity.jpg';
outputFilename3 = 'g2.jpg';

nMax = 10; % Maximum matrix dimension
nVector = (-1:nMax+1)';
rho = zeros(nMax+3,1);
rho(nMax+2) = 1;
h = 0.1;
t = 0:h:100;
I = zeros(length(t),1);
g2 = zeros(length(t),1);
rhoSaved = zeros(nMax+3,length(t));


% Solve with Runge-Kutta
for j = 1:length(t)
    I(j) = nVector' * rho;
    g2(j) = (nVector.*(nVector-1))'*rho/I(j)^2;
    rhoSaved(:,j)=rho;
    
    rho = RungeKutta(h, rho, @rhoPoint );
    
    rho = rho / sum(rho);
end

for i = 1:nMax+2
   plot(t,rhoSaved(i,:),'lineWidth',1.5,'DisplayName', ...
       strcat('$\rho_{',num2str(i-1),'}$'));
   hold on;
end
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultTextInterpreter','latex');
legend ('Location','northeast');
xlabel('t [s]');
ylabel('Probability [a.u.]');
print(outputFilename,outputFiletype);

hold off;
plot(t,I);
xlabel('t [s]');
ylabel('Intensity [a.u.]');
print(outputFilename2,outputFiletype);

hold off;
plot(t,g2);
xlabel('t [s]');
ylabel('g^{(2)} [a.u.]');
print(outputFilename3,outputFiletype);
end

function deriv = rhoPoint( rho)
    Gamma = 0.1; % Coupling constant
    nB = 1; % 1/(exp(beta*Omega)-1) Bose-Einstein bath occupation
    nMax = 10; % Maximum matrix dimension
    nVector = (-1:nMax+1)';
    deriv = zeros(length(rho),1);
    deriv(2:end-1) = Gamma .* (...
        nB * nVector(2:end-1) .* rho(1:end-2) - ...
        (nVector(2:end-1)+(2*nVector(2:end-1)+1)*nB).*rho(2:end-1) + ...
        (1+nB)*nVector(3:end).*rho(3:end) );
    
end