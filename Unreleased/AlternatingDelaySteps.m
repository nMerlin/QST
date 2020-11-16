Start = -115;
End = 175;
Increment = 0.8;

N = floor((End - Start)/Increment/2);
k = 0:N;
x1 = Start + k*Increment;
x2 = End - k*Increment;

x = cat(1,x1,x2);
x = x(:)';