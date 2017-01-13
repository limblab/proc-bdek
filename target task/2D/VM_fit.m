function[pfits,thetas,curve,base,gain,fwhm] = VM_fit(x,y,theta_list,varargin)
% function[pfits,thetas,curve,base,gain,fwhm] = VM_fit(x,y,theta_list,varargin)

xs = reshape(x,[],1);
ys = reshape(y,[],1);

%Initial cosine fit estimates
b = [ones(size(xs)) cos(xs) sin(xs)]\ys;

pd_i = mod(atan2(b(3),b(2))+10*pi,2*pi);
gain_i = sqrt(b(2).^2+b(3).^2);
off_i = b(1)-gain_i;
k_i = pi/4;

b_i = [off_i gain_i k_i pd_i];
% bounds_L = [0 -100 0.1 -2*pi];
% bounds_H = [100 100 25 4*pi];
bounds_L = [0 -100 0.1 -2*pi];
bounds_H = [100 100 25 4*pi];

vm_func = @(xs,p) p(1) + p(2)*exp(p(3)*cos(xs-p(4)));

costfunc = @(p) sum((ys - vm_func(xs,p)).^2);

fmopts = optimset('Algorithm','interior-point','Display','off','MaxIter',10000);
%[pfits,~,exitflag] = fminsearch(costfunc,b_i,fmopts);
[pfits,~,exitflag] = fmincon(costfunc,b_i,[],[],[],[],bounds_L,bounds_H,'nonlin_helper',fmopts);

if nargin > 2
    thetas = theta_list;
else
    thetas = 0.01:0.01:2*pi;
end

if exitflag==0
    fprintf('no convergence\n');
end
    
curve = vm_func(thetas,pfits);

base = min(curve);
gain = max(curve)-min(curve);

split = curve - 0.5*(min(curve)+max(curve));
fwhm = sum(split >= 0).*diff(thetas(1:2));

% outs = [curve base gain fwhm];

end



