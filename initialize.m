function [x0,A,b,actcon] = initialize(model)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[A,b,xlb,xub] = constraints(model);
x0 = initguess(model,A,b,xlb,xub);
A = -A;
b = -b;
nx = length(x0);
A = [A;eye(nx);-eye(nx)];
b = [b;xlb;-xub];
actcon = find(A*x0<=b);



end

