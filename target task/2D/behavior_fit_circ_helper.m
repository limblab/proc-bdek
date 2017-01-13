function [k1_k2_k3] = behavior_fit_circ_helper(mu1_mu2_reaches)

% Functions
mufunc = @(mu1,mu2,kratio) mu1 + atan(sin(mu2-mu1)./(kratio + cos(mu2-mu1)));
mufitfunc = @(mu1,mu2,reaches,kratio) sum(circ_dist(mu1 + atan(sin(mu2-mu1)./(kratio + cos(mu2-mu1))),reaches).^2);
k2func = @(mu1,mu2,k3,kratio) sqrt(k3^2./(kratio^2 + 2*kratio*cos(mu2-mu1) + 1));

% Assign
mu1 = mu1_mu2_reaches(:,1);
mu2 = mu1_mu2_reaches(:,2);
reaches = mu1_mu2_reaches(:,3);

kratfit = fminbnd(@(x) mufitfunc(mu1,mu2,reaches,x),eps,1000);

estreaches = mufunc(mu1,mu2,kratfit);

k3 = circ_kappa(circ_dist(estreaches,reaches));
k2 = mean(k2func(mu1,mu2,k3,kratfit));
k1 = kratfit*k2;

k1_k2_k3 = [k1 k2 k3];
k1_k2_k3(k1_k2_k3<0.5)=.5;
end