function [c,ceq] = nonlin_helper(x)

c =  -(x(1) + x(2)*exp(-x(3)));
ceq = [];

end