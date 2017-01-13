function[k_est, k_low, k_high] = kappa_calc(samps,a)

%Based on "Improved efficient approximation of concentration parameter and
%confidence interval for circular distribution"
% Hassan, Hussin, Zubairi: ScienceAsia (2012)
n = length(samps);
C = 1/n*sum(cos(samps));
S = 1/n*sum(sin(samps));

R = (C^2 + S^2);
t = sqrt(R);

if t < 0.6137

    sols = roots([1 0 -6 0 48 -96*t]);
    k_est = sols(imag(sols)==0);

elseif t >= 0.6137

    sols = roots([(8*t - 8) 4 1 1]);
    k_est = sols(imag(sols)==0);

else
    fprintf('Yo'' samples are messed up, son!\n');
end

%% Confidence %%%
% Method 1

v = sqrt(-2*log(t));

Y = exp(-(n-1)*v^2/(2*chi2inv(1-a/2,n-1)));
Z = exp(-(n-1)*v^2/(2*chi2inv(a/2,n-1)));

if Y < 0.6137

    sols = roots([1 0 -6 0 48 -96*Y]);
    l_k = sols(imag(sols)==0);

elseif Y >= 0.6137

    sols = roots([(8*Y - 8) 4 1 1]);
    l_k = sols(imag(sols)==0);

else
    fprintf('Yo'' samples are messed up, son!\n');
end

if Z < 0.6137

    sols = roots([1 0 -6 0 48 -96*Z]);
    h_k = sols(imag(sols)==0);

elseif Z >= 0.6137

    sols = roots([(8*Z - 8) 4 1 1]);
    h_k = sols(imag(sols)==0);

else
    fprintf('Yo'' samples are messed up, son!\n');
end

k_low1 = h_k;
k_high1 = l_k;

%% Method 2

B = 1.96/sqrt(n*(1 - t/k_est - t^2));

k_low2 = -B + k_est;
k_high2 = B + k_est;

%% Output
if (nargin < 3 || (nargin > 2 && method == 1))
    k_low = k_low1;
    k_high = k_high1;
elseif nargin > 2 && method == 2
    k_low = k_low2;
    k_high = k_high2;
end
   







