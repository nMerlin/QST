function [Pn]= photonNumberDist(X, filename)
%  Reconstructs the Photon Number Distribution (the diagonal part of the 
%  density matrix) from a phase-averaged state.
%  It only delivers proper results for mean photon numbers up to ca 50
%  due to calculation of Hermite Polynoms.
%   
%   Input Parameters:
%   X - Matrix containing the measured quadrature values. You can create it
%       for example with the function PREPAREPHAVDATA
%   FILENAME - this string will be included in the filename of the png-file


%Remove offset
X = X(:);
X = X - mean(mean(X));

% Mean photon number
nAv = mean(X.^2)-0.5;

uniq = unique(X(:));
maxValue = max(-min(uniq),max(uniq));
hDisc = min(diff(uniq)); % discretization
histEdges = (-maxValue-hDisc/2):hDisc:(maxValue+hDisc/2);
[Pq,~]=histcounts(X,histEdges,'Normalization','probability');
q = -maxValue:hDisc:maxValue; 

%max. photonnumber
Nmax = min(max(ceil(4*nAv),10),95);

%detection efficiency
eff = 1;

%Calculate Hermite-matrix
H = zeros(length(q),Nmax+1);
for n = 0:Nmax
    H(:,n+1) = polyval(hermitePoly(n),q) .* exp(-q.^2).* polyval(hermitePoly(n),q) ;
    %H(:,n+1) = polyval(hermitePoly(n),q).^2 .* exp(-q.^2);
    %H(:,n+1) = 2*log(polyval(hermitePoly(n),q)) -q.^2;  %with log
end

if eff == 1
    %Calculate Coefficients vector
    a = zeros(Nmax+1,1);
    a(1) = 1/sqrt(pi);
    for n = 1:Nmax
        a(n+1) = a(n)/(2*(n));
    end

    %Calculate Transfer-matrix
    A = zeros(length(q),Nmax+1);
    for n = 1:Nmax+1
        A(:,n) = H(:,n)*a(n);
        %A(:,n) = H(:,n) + log(a(n)); %with log
    end
     
else

    %Calculate Coefficients Matrix  
    Acoeff = zeros(Nmax+1,Nmax+1);
    for n = 0:Nmax
        for m = 0:n
            Acoeff(n+1,m+1) = (1-eff)^(n-m)*eff^m*factorial(n)/(sqrt(pi)*2^m*factorial(n-m)*(factorial(m))^2);  
        end
    end

    %Calculate Transfer-matrix
    A = zeros(length(q),Nmax+1);
     for n = 1:Nmax+1
         for m = 1:n
             A(:,n) = A(:,n) + H(:,m).*Acoeff(n,m);
         end     
     end
end

%A = exp(A);

Pn = A\Pq';   

Pn = Pn./sum(Pn);

n = 0:Nmax;
nAv2 = n*Pn;

bar(n,Pn);
xlabel('n','fontsize',13);
ylabel('\rho_{nn}','fontsize',13);
title({filename,['Phase-averaged quantum state (n=',num2str(nAv),')']}, ...
    'HorizontalAlignment','center', 'FontWeight','bold');
text('Units','normalized','Position',[0.97 0.95],'String',...
    ['n (\rho_{nn})= ',num2str(nAv2)], 'HorizontalAlignment','right','Interpreter', 'tex');
print(strcat('nDist-',filename,'-',num2str(nAv),'photons.png'), '-dpng');


end