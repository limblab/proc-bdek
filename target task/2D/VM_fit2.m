function[pfits,thetas,curve,base,gain,fwhm] = VM_fit2(x,y,theta_list,varargin)
% function[pfits,thetas,curve,base,gain,fwhm] = VM_fit(x,y,theta_list,varargin)

xs = reshape(x,[],1);
ys = reshape(y,[],1);

%Initial cosine fit estimates
b = [ones(size(xs)) cos(xs) sin(xs)]\ys;

pd_i = mod(atan2(b(3),b(2))+10*pi,2*pi);
gain_i = sqrt(b(2).^2+b(3).^2);
off_i = b(1)-gain_i;
k_i = 1;

if off_i>0; binit1 = off_i; else binit1 = 0; end
if (off_i+gain_i)>0; binit2 = (off_i+gain_i); else binit2 = 0; end

b_i = [binit1 gain_i k_i];% pd_i];
b_i2 = [binit2 , -gain_i, k_i];% pd_i+pi];
% bounds_L = [0 -100 0.1 -2*pi];
% bounds_H = [100 100 25 4*pi];
bounds_L = [0   0   0.01];% -2*pi];
bounds_H = [100 100 25];% 4*pi];

bounds_L2 = [0  -100 0.01];% -2*pi];
bounds_H2 = [100 0   25];%   4*pi];

% vm_func = @(xs,p) p(1) + p(2)*exp(p(3)*cos(xs-p(4)));

vm_func1 = @(xs,p) p(1) + p(2)*exp(p(3)*cos(xs));
vm_func2 = @(xs,p) p(1) + p(2)*exp(p(3)*cos(xs-pi));

% costfunc = @(p) sum((ys - vm_func(xs,p)).^2);

costfunc1 = @(p) sum((ys - vm_func1(xs,p)).^2);
costfunc2 = @(p) sum((ys - vm_func2(xs,p)).^2);


% fmopts = optimset('Algorithm','interior-point','Display','off','MaxIter',10000);
fmopts = optimset('Algorithm','active-set','Display','off','MaxIter',10000);
%[pfits,~,exitflag] = fminsearch(costfunc,b_i,fmopts);
[pfits1,~,exitflag] = fmincon(costfunc1,b_i,[],[],[],[],bounds_L,bounds_H,@(x)nonlin_helper(x),fmopts);

ypreds = reshape(vm_func1(x,pfits1),[],1);
SSE_1 = sum((ypreds-ys).^2);

[pfits2,~,exitflag2] = fmincon(costfunc2,b_i2,[],[],[],[],bounds_L2,bounds_H2,@(x)nonlin_helper(x),fmopts);

ypreds2 = reshape(vm_func2(x,pfits2),[],1);
SSE_2 = sum((ypreds2-ys).^2);

if nargin > 2
    thetas = theta_list;
else
    thetas = 0.01:0.01:2*pi;
end
 

if SSE_1 < SSE_2
    pfits = pfits1;
    %if exitflag==0; fprintf('no convergence\n'); end 
    curve = vm_func1(thetas,pfits1);
else
    pfits = pfits2;
    %if exitflag2==0; fprintf('no convergence\n'); end 
    curve = vm_func2(thetas,pfits2);
end

%curve = vm_func(thetas,pfits);

base = min(curve);
gain = max(curve)-min(curve);

split = curve - 0.5*(min(curve)+max(curve));
fwhm = sum(split >= 0).*diff(thetas(1:2));

% outs = [curve base gain fwhm];

end



