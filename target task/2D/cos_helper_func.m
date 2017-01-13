function[curve] = cos_helper_func(xs,ys)

b = glmfit([cos(xs) sin(xs)],ys,'poisson');

thets = -pi:0.01:pi;

curve = glmval(b,[cos(thets') sin(thets')],'log');