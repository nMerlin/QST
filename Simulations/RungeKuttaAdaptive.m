function [ynew] = RungeKuttaAdaptive(h, yold, deriv)
%First version of an adaptive RungeKutta Solver. Input is stepsize h, vector with old y
%values, and function handle to compute the derivation from yold. Output ist new y. 

yscal = yold; %Might be chosen differently
eps = 10^-6;  
delta0 = eps * yscal; %desired error

y1 = RK(0.5 * h, yold, deriv);  %half step
y2 = RK(h, yold, deriv); %full step

delta = abs(y1 - y2); %actual error

if mean(delta) > mean(delta0)    %actual error too large
    h = h * (mean(delta0)/mean(delta))^0.2;   %compute better stepsize
    %disp(h);
end

ynew = RK(h, yold, deriv);

end

function [ynew] = RK(h, yold, deriv)

k1 = h * deriv(yold);
k2 = h * deriv(yold + 0.5 * k1);
k3 = h * deriv(yold + 0.5 * k2);
k4 = h * deriv(yold + k3);
ynew = yold + (k1 + 2*k2 + 2*k3 + k4)/6;


end