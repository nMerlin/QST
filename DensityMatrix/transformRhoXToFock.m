function[rho] = transformRhoXToFock(rhoX,nMax,xy)
% transforms a density matrix rhoX given in spatial coordinate base to the
% Fock base with maximum photon number nMax.
% xy is the coordinate axis corresponding to rhoX, e.g. -5:0.125:5
% see Lvovsky 2009 or d'Ariano 1994, eq.3

nX = length(xy);

%% make x and y matrix
xMatrix = zeros(nX,nX);
for index = 1:nX
    xMatrix(index,:) = xy;
end
yMatrix = xMatrix';

rho = zeros(nMax+1,nMax+1);
for n = 0:nMax   
    for m = 0:nMax
        Conv = rhoX.*fockstate(n,xMatrix).*fockstate(m,yMatrix);
        rho(n+1,m+1) = sum(sum(Conv));        
    end;
end;
rho = normalize(rho); % normalize to unity trace

end