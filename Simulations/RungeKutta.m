function [yold] = RungeKutta(h, yold, deriv)

k1 = h * deriv(yold);
k2 = h * deriv(yold + 0.5 * k1);
k3 = h * deriv(yold + 0.5 * k2);
k4 = h * deriv(yold + k3);
ynew = yold + (k1 + 2*k2 + 2*k3 + k4)/6;


end
