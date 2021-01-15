function [a,b,c] = doubleToUDP(x)

a = sign(x) +1;
b = floor(abs(x));
c = (round((abs(x) - b),2))*100;
