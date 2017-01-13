function[PD,b] = PD_helper_func_COS(xs_ys_spds)

xsys = xs_ys_spds(:,1:2);

xs = xsys(:,1);
ys = xsys(:,2);

A = [ones(size(xs)), cos(xs), sin(xs)]; 
b = A\ys;

PD = atan2(b(3),b(2));
end