function [expQ,expP] = postselection(pps,qps,n,transmission)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

t = sqrt(transmission);
r = sqrt(1-transmission);
expQ = (n*qps*r*t)/(pi*(-1+n*(-1+t^2))^2)*exp((pps^2+qps^2)/(-1+n*(-1+t^2)));
expP = expQ;

end

