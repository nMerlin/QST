function [ I ] = cavityThermalBath_lin()
%Simulates a single mode in a cavity with a surrounding bath
%Computes time evolution as an eigenvalue problem

% Parameters
outputFilename = 'CavityThermalBath_lin.jpg';
outputFiletype = '-djpeg';
outputFilename2 = 'Intensity_lin.jpg';
outputFilename3 = 'g2_lin.jpg';

Gamma = 0.1; % Coupling constant
nB = 1; % 1/(exp(beta*Omega)-1) Bose-Einstein bath occupation

nMax = 10; % Maximum photon number
nVector = (0:nMax)';
rho = zeros(nMax+1,1);
rho(nMax+1) = 1;

%Compute Transfer Matrix

diagmain = -Gamma*(nVector + (2*nVector + 1)*nB);
diagleft = Gamma * nB * nVector(2:end);
diagright = Gamma* (1 + nB) * (nVector(1:end-1)+1);

M = diag(diagmain,0) + diag(diagleft,-1) + diag(diagright,1);

h = 0.1;
t = 0:h:100;
I = zeros(length(t),1);
g2 = zeros(length(t),1);

rhoSaved = zeros(nMax+1,length(t));

%Compute Eigenvalues and Eigenvectors
[V,D] = eig(M); %Eigenvectors V, diagonal matrix with eigenvalues D
%decomposite starting rho into eigenvectors
b = V\rho; %vector b with coefficients

%Compute time evolution
for j = 1:length(t)
    
    I(j) = nVector' * rho;
    g2(j) = (nVector.*(nVector-1))'*rho/I(j)^2;
 
    rhoSaved(:,j)=rho;
    
    ex = exp( D * t(j)*ones(nMax+1,1));
    rho = V * (b.*ex);
    rho = rho / sum(rho);

end

%compute assumed stationary state
a = nB/ (1+nB); % a is exp(-beta*omega)
rho0 = zeros(nMax + 1,1);
rho0(1) = (1-a)/(1-a^(nMax+1));
for i = 2:nMax + 1
    rho0(i) = rho0(i-1)*a;
end


for i = 1:nMax+1 
   plot(t,rhoSaved(i,:),'lineWidth',1.5,'DisplayName', ...
       strcat('$\rho_{',num2str(i-1),'}$'));
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
axis tight
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




