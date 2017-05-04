function [ yout, xout ] = testStepperDopr853( )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

htry = 0.01;
yout = zeros(length(0:htry:1),1);
xout = yout;
y = 1;
x = 0;
atol = 1e-6;
rtol = 1e-6;
dense = false;
Stepper = StepperDopr853(y, derivs(x,y), x, atol, rtol, dense);
Stepper.step(htry, @derivs);

for iYout = 1:length(yout)
    yout(iYout) = Stepper.y;
    xout(iYout) = Stepper.x;
    %x = (j-1)*h;
    Stepper.step(Stepper.hnext, @derivs);
    %y = tooSimpleRK4(y, derivs(x,y), x, h, @derivs);
end

% Compare with exact solution
x = 0:0.01:Stepper.x;
y = exp(-x.^2);
plot(xout,yout);
hold on;
plot(x,y);

end

function dydx = derivs(x,y)
    dydx = -2*x*y;
end