function[Bayes_est, posterior, post_k] = Bayes_circ_trial(mu_likelihood,k_likelihood,mu_prior,k_prior)

% DEFINE FUNCTIONS TO MULTIPLY TWO VON MISES
% From "Cue combination on the circle and the sphere" (Richard F. Murray
% and Yaniv Morgenstern)
get_mu3_k3 = @(mu1,mu2,k1,k2) [mu1 + atan2(sin(mu2-mu1),(k1/k2)+cos(mu2-mu1)), ...
                               sqrt(k1^2 + k2^2 + 2*k1*k2*cos(mu2-mu1))];                       
prodvm = @(k1,k2,theta,mu3_k3) besseli(0,mu3_k3(2))/(2*pi*besseli(0,k1)*...
                               besseli(0,k2))*circ_vmpdf(theta,mu3_k3(1),mu3_k3(2));
                           
% variables                     
thetas = 0:0.01:2*pi;

%%%%%%%%%%%%%%%%%%% POSTERIOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiply likelihood and prior
m3k3 = get_mu3_k3(mu_likelihood,mu_prior,k_likelihood,k_prior); 
posterior = prodvm(k_likelihood,k_prior,thetas,m3k3);

Bayes_est = m3k3(1);
post_k = m3k3(2);

