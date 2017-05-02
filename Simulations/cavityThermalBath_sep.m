function [ I ] = cavityThermalBath_sep()
%Simulates a single mode in a cavity with a surrounding bath

% Parameters
outputFilename = 'CavityThermalBath.jpg';
outputFiletype = '-djpeg';
outputFilename2 = 'Intensity.jpg';
outputFilename3 = 'g2.jpg';
outputFilename4 = 'h.jpg';

nMax = 10; % Maximum matrix dimension
nVector = (-1:nMax+1)';
rho = zeros(nMax+3,1);
rho(nMax+2) = 1;
h = 0.1;
t = 0:h:100;
I = zeros(length(t),1);
g2 = zeros(length(t),1);
hsaved = zeros(length(t),1);
rhoSaved = zeros(nMax+3,length(t));
err = 10^-6;

% Solve with Runge-Kutta
for j = 1:length(t)
    I(j) = nVector' * rho;
    g2(j) = (nVector.*(nVector-1))'*rho/I(j)^2;
    hsaved(j) = h;
    rhoSaved(:,j)=rho;
    
    %rho = RungeKutta(h, rho, @rhoPoint );  %no adaptive stepsize
    %[rho,hnew] = RungeKuttaAdaptive(h, rho, @rhoPoint );  %my Adaptive
    [rho, h] = rka(rho,h,err,@rhoPoint);   %from matlab file exchange
    rho = rho / sum(rho);
    hsaved(j) = h;
end


%compute assumed stationary state
nB = 1;
a = nB/ (1+nB); % a is exp(-beta*omega)
rho0 = zeros(nMax,1);
rho0(1) = (1-a)/(1-a^(nMax+1));
for i = 2:nMax
    rho0(i) = rho0(i-1)*a;
end
disp(rho0);

%Plotting

for i = 2:nMax+2
   plot(t,rhoSaved(i,:),'lineWidth',1.5,'DisplayName', ...
       strcat('$\rho_{',num2str(i-2),'}$'));
   hold on;
end


plot(t(end)*ones(length(rho0),1),rho0,'o','MarkerFaceColor', 'k', 'MarkerEdgeColor','k', ...
    'LineWidth', 2, 'MarkerSize', 4,'DisplayName', ...
       strcat('$\rho_{stat}$'));
hold on;

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

hold off;
plot(t,hsaved);
xlabel('t [s]');
ylabel('h [a.u.]');
print(outputFilename4,outputFiletype);
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





