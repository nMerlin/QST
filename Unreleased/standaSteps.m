xmin = 163498;
xmax=166798;
number = 100;
Prange = 100;
del = 2*(xmax-xmin);
k = 0:number-2;
a = Prange^(1/length(k));  
x = xmin + del/pi * asin(sqrt(1/Prange * a.^k));
x = round(x);
x = [x xmax];
P = (sin(pi/del*(x-xmin))).^2;

%% to have more data points for high powers 
xmin = 163400;
xmax=166798;
number = 75;
Prange = 100;
del = 2*(xmax-xmin);
k = [0:round(number/3)  round(number/3):0.5:round(number*2/3 -2)];
a = Prange^(1/(number*2/3-2));  
x = xmin + del/pi * asin(sqrt(1/Prange * a.^k));
x = round(x);
x = [xmin x xmax];
P = (sin(pi/del*(x-xmin))).^2;