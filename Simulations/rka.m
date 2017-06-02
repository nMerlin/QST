function [xSmall, tau] = rka(x,tau,err,func)
% Adaptive Runge-Kutta routine from https://de.mathworks.com/matlabcentral/fileexchange/55830-adaptive-runge-kutta-integrator?focused=6221084&tab=function
% Inputs
% x Current value of the dependent variable
% t Independent variable (usually time)  ; deleted because not needed 
% tau Step size (usually time step)
% err Desired fractional local truncation error
% func Right hand side of the ODE; func is the
% name of the function which returns dx/dt
% Calling format func(x,t).
% Outputs
% xSmall New value of the dependent variable
% t New value of the independent variable
% tau Suggested step size for next call to rka
%* Set initial variables
xSave = x;
% Save initial values
safe1 = .9;  %.9
safe2 = 4.;
% Safety factors
%* Loop over maximum number of attempts to satisfy error bound
maxTry = 100;
for iTry=1:maxTry
    %* Take the two small time steps
    half_tau = 0.5 * tau;    
    xTemp = rk4(xSave,half_tau,func);
    xSmall = rk4(xTemp,half_tau,func);
    %* Take the single big time step
    xBig = rk4(xSave,tau,func);
    %* Compute the estimated truncation error
    scale = err * (abs(xSmall) + abs(xBig))/2.;
    xDiff = xSmall - xBig;
    errorRatio = max( abs(xDiff)./(scale + eps) );
    %* Estimate new tau value (including safety factors)
    tau_old = tau;
    tau = safe1*tau_old*errorRatio^(-0.20);
    tau = max(tau,tau_old/safe2);
    tau = min(tau,safe2*tau_old);
    %* If error is acceptable, return computed values
    if (errorRatio < 1)
        return
    end
end
%* Issue error message if error bound never satisfied
% error('ERROR: Adaptive Runge-Kutta routine failed');
% return;

end

function [ynew] = rk4(yold, h, deriv)

k1 = h * deriv(yold);
k2 = h * deriv(yold + 0.5 * k1);
k3 = h * deriv(yold + 0.5 * k2);
k4 = h * deriv(yold + k3);
ynew = yold + (k1 + 2*k2 + 2*k3 + k4)/6;


end