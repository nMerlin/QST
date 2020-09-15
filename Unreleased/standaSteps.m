x0 = 163498;
xmax=166798;
del = 2*(xmax-x0);
k = 0:29;
Prange = 1000;
a = Prange^(1/length(k));
x = x0 + del/pi * asin(sqrt(1/Prange * a.^k));
x = round(x);
x = [x xmax];

P = (sin(pi/del*(x-x0))).^2;