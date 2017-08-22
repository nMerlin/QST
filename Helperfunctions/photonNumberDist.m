function [Pn]= photonNumberDist(X, filename)
%  Reconstructs the Photon Number Distribution (the diagonal part of the 
%  density matrix) from a phase-averaged state.
%  It only delivers proper results for mean photon numbers up to ca 50
%  due to calculation of Hermite Polynomials.
%  The idea comes from 
%  Konrad Banaszek: Maximum likelihood estimation of photon number distribution
%  from homodyne statistics, Phys. Rev. A 57, 5013-5015 (1998)
%     (arXiv:physics/9712043)
%  According to this, the phase-averaged probability distribution p(q) 
%  is a linear combination of the photon number distribution p_n:
%  P(q) = sum_{n=0}^{\infty} A_n (q) * p_n  (1), 
%  with A_n (q) = sum_{m=0}^{n} (1-eff)^(n-m) * eff^m * n! /(sqrt(pi)* 2^m * (n-m)! * (m!)^2) 
%  * (H_m (q))^2 * exp(-q^2).  (2)
%  Eff is the detection efficiency, H_m is the mth Hermite Polynomial. 
%  The code cuts the sum in (1) off at n = Nmax and solves the linear equation system for p_n.
%  In the paper, a maximum likelihood method is proposed. 

%   Input Parameters:
%   X - Matrix containing the measured quadrature values. You can create it
%       for example with the function PREPAREPHAVDATA
%   FILENAME - this string will be included in the filename of the png-file
%   The function hermitePoly(n) is needed.

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
    end
     
else

    %Calculate Coefficients Matrix  
    Acoeff = zeros(Nmax+1,Nmax+1);
    for n = 0:Nmax
        for m = 0:n
            Acoeff(n+1,m+1) = (1-eff)^(n-m)*eff^m*factorial(n)/ ...
                (sqrt(pi)*2^m*factorial(n-m)*(factorial(m))^2);  
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