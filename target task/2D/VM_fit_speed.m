function[pfits,thetas,curve,base,gain,fwhm] = VM_fit_speed(x,y,speed,theta_list,varargin)
% function[pfits,thetas,curve,base,gain,fwhm] = VM_fit3(x,y,theta_list,varargin)

xs = reshape(x,[],1);
ys = reshape(y,[],1);
sps = reshape(speed,[],1);
tlist = reshape(theta_list,[],1);

[b,devfull] = glmfit([sps.*cos(xs) sps.*sin(xs)],ys,'poisson');
[~,devnull] = glmfit(ones(size(ys)),ys,'poisson','constant','off');

pval = 1 - chi2cdf(devnull-devfull,2);

%Initial cosine fit estimates
if nargin > 3
    thetas = theta_list;
else
    thetas = 0.01:0.01:2*pi;
end
  
if pval < 0.1 && sqrt(b(2).^2+b(3).^2)<1e3
    curve = glmval(b,[cos(tlist) sin(tlist)],'log')';
    pfits = b;
else
    curve = nan(1,length(tlist));
    pfits = NaN(size(b));
end

base = min(curve);
gain = max(curve)-min(curve);

split = curve - 0.5*(min(curve)+max(curve));
fwhm = sum(split >= 0).*diff(thetas(1:2));
end
%%