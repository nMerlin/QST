function [pattern] = patternFunction(XGrid, R)
 % compute the pattern function from the X parameters, see Jans notes, equ. 14 and 15
    fun = @(u,XGridf,Rf) u .* (acos(u) - u.*sqrt(1 - u.^2)) .* cos(u.*XGridf).*exp(2*Rf^2 *u.^2);
    h = integral(@(u)fun(u,XGrid,R),0,1,'ArrayValued',true );   
    pattern = (16*R^2 / pi^3).*h;


end