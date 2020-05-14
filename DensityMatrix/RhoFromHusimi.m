function [rho,rhoX] = RhoFromHusimi(HS,varargin)
% DOES NOT WORK CORRECTLY!!!
% computes the density matrix from the Husimi Q function (HS). HS must
% be a quadratic matrix. 
% Formula from equation (7) and (8) from V. A. Andreev: A TRANSFORMATIONAL
% PROPERTY..., Theoretical and Mathematical Physics, 166(3): 356–368 (2011)
% https://link.springer.com/content/pdf/10.1007/s11232-011-0028-8.pdf
% Optional Input Arguments:
% 'PQ': Specify q axis, which is also x axis. Default is -20:0.125:20, 
% as this is used in mainWigner. It must have the same length as the HS
% matrix.
% 'nMax': maximum photon number for fock base.
% 'Method': which equation is used, can be 'eq7', 'eq8-1' or 'eq8-2'.
% 'Filename' for saved data

%% Validate and parse input arguments
parser = inputParser;
defaultPQ = -20:0.125:20;
addParameter(parser,'PQ',defaultPQ,@isvector);
defaultNorm = 1/sqrt(2); % sqrt([q,p]/(2*i))
addParameter(parser,'Norm',defaultNorm,@isnumeric);
defaultNMax = 30; 
addParameter(parser,'nMax',defaultNMax,@isnumeric);
defaultMethod = 'eq8-1'; 
addParameter(parser,'Method',defaultMethod);
defaultFilename = ''; % sqrt([q,p]/(2*i))
addParameter(parser,'Filename',defaultFilename);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[filename,method,nMax,norm,pq] = c{:};


hbar = 1;
%hbar = 1.054571817e-34;
%% prepare p, q axis
[p,q, x, y] = deal(pq);
nP = length(p);
nHS = length(HS);
shift = (nHS-nP)/2;
HSs = HS(shift+1:end-shift,shift+1:end-shift);

%% make p and q matrix
qMatrix = zeros(nP,nP);
for index = 1:nP
    qMatrix(index,:) = q;
end
pMatrix = qMatrix';

%% prepare Density matrix in coordinate space. The coordinates will also use the pq axis. 
rhoX = zeros(size(HSs));
[maxX,maxY]=size(HSs);

%% evaluate each value of the Density matrix in coordinate space.
dispstat('start rhoX','timestamp','keepthis',0);
tic;
% Method from Andreev 2001, equ. 7; 
if strcmp(method,'eq7')
    for xi = 1:maxX   
        for yi = 1:maxY
            f1 = exp((x(xi)-y(yi))^2/2);
            Kernel = 1/sqrt(pi) * exp(-(qMatrix - 0.5*(x(xi)+y(yi))).^2 - (x(xi)-y(yi))^2 + 1i*pMatrix*(x(xi)-y(yi)));
            Prod = Kernel .* HSs;
            Sum = 0;
            for n = 0:nMax % not sure if it is ok to stop the sum at nMax...
                     Conv = polyval(hermitePoly(2*n),q-0.5*(x(xi)+y(yi))^2) * Prod;
                     Sum = Sum + (-1)^n/(2^n*factorial(n))  *  sum(Conv);

    %             Conv = polyval(hermitePoly(2*n),qMatrix-0.5*(x(xi)+y(yi))^2) .* Kernel .* HSs;  
    %             Sum = Sum + (-1)^n/(2^n*factorial(n))  *  sum(sum(Conv)); 
            end
            rhoX(xi,yi) = f1 * Sum;        
        end;
    end;

% Method from Andreev 2001, equ. 8)
elseif strcmp(method,'eq8-1')
    for xi = 1:maxX   
        for yi = 1:maxY
            f1 = 1i/sqrt(pi)*exp((x(xi)^2+y(yi)^2)/2);
            Kernel = exp(hbar^2*qMatrix.^2 +1i*pMatrix*(-x(xi) + y(yi))  + hbar*qMatrix*(-x(xi) - y(yi)));        
            Prod = Kernel .* HSs;
            rhoX(xi,yi) = f1 * sum(sum(Prod));        
        end;
    end;

% Method from Andreev 2001, equ. 8), different coordinates; gives the
% same result as the above method 
elseif strcmp(method,'eq8-2') 
    p1Matrix = pMatrix - 1i*hbar*qMatrix;
    p2Matrix = pMatrix + 1i*hbar*qMatrix;
    for xi = 1:maxX   
        for yi = 1:maxY
            f1 = 1/(2*hbar)*1/sqrt(pi)*exp((x(xi)^2+y(yi)^2)/2);
            Kernel = exp(-1/4 *(p1Matrix-p2Matrix).^2 -1i*p1Matrix*x(xi) + 1i*p2Matrix*y(yi));        
            Prod = Kernel .* HSs;
            rhoX(xi,yi) = f1 * sum(sum(Prod));        
        end;
    end;
end;
toc;
dispstat('end rhoX','timestamp','keepthis',0);
rhoX(isnan(rhoX))=0; % remove nan values
rhoX = normalize(rhoX); % normalize to unity trace

%% transfer density matrix to the Fock base.
dispstat('start rhoFock','timestamp','keepthis',0);

rho = transformRhoXToFock(rhoX,nMax,pq);

dispstat('end rhoFock','timestamp','keepthis',0);
save(['rho-', method, filename,'.mat'],'rho','rhoX');

end