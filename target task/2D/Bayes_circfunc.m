function[Bayes_est, posterior, post_k,mu_likelihood] = Bayes_circfunc(tt,slices)
% Bayes_circ(tt,slices,priortype,trialwindow,varargin) returns the optimal
% Bayesian estimate (Bayes_est) as well as the full posterior distribution
% for each trial (posterior) and the concentration metric of the posterior 
% (post_k). 
%
% tt: trial table in which:
%     column  2: target location
%
% slices: (n x c) matrix containing locations of c target cues for each 
%         trial n
%
% priortype: specifies the type of prior to use
%          - 'all'     : uses the mean and variance from all trials
%          - 'running' : for each trial uses the previously seen target
%                        locations
%          - 'history' : Uses only the past 'trialwindow' number of trials
%                        to compute prior. 
%
% trialwindow: see priortype: 'history'    


% DEFINE FUNCTIONS TO MULTIPLY TWO VON MISES
% From "Cue combination on the circle and the sphere" (Richard F. Murray
% and Yaniv Morgenstern)
get_mu3_k3 = @(mu1,mu2,k1,k2) [mu1 + atan2(sin(mu2-mu1),(k1/k2)+cos(mu2-mu1)), ...
                               sqrt(k1^2 + k2^2 + 2*k1*k2*cos(mu2-mu1))];          
prodvm = @(k1,k2,theta,mu3_k3) besseli(0,mu3_k3(2))/(2*pi*besseli(0,k1)*...
                               besseli(0,k2))*circ_vmpdf(theta,mu3_k3(1),mu3_k3(2));
                           
% variables                     
thetas = 0:0.01:2*pi;
nt = size(tt,1);

% initialize
posterior = NaN(nt,length(thetas));
[Bayes_est, post_k] = deal(NaN(nt,1));
pflag = zeros(nt,1);
mu_likelihood = zeros(size(tt,1),1);

k_prior = circ_kappa(tt(:,2)); 
mu_prior = circ_mean(tt(:,2));
for i = 1:nt

    %%%%%%%%%%%%%%%%%%% LIKELIHOOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lines = slices(i,:);
    ns = length(lines);
    
    k_likelihood = circ_kappa(lines)*ns;
    mu_likelihood(i) = circ_mean(lines');
%     
    %%%%%%%%%%%%%%%%%%% PRIOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     if i > 1
%         
%         
%         %KP = k_prior + G*tan(h./2);
%         %KP = k_prior + G*(2/(1+exp(-h)) - 1);
%         KP = k_prior + G(likes==tt(i,3))*h;
%     else
%         KP = k_prior;
%     end
    KP = k_prior;

    %%%%%%%%%%%%%%%%%%% POSTERIOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Multiply likelihood and prior
    m3k3 = get_mu3_k3(mu_likelihood(i),mu_prior,k_likelihood,KP); 
    posterior(i,:) = prodvm(k_likelihood,KP,thetas,m3k3);
    post_k(i) = m3k3(2);

    if isnan(m3k3(1))
        Bayes_est(i) = mu_likelihood(i); 
        pflag(i) = 1;
    else
        Bayes_est(i) = m3k3(1);
    end

end

%Output warnings if necessary
if sum(pflag > 0);   
    fprintf('Ill-defined estimate for trials: ');
    fprintf('%d ',find(pflag==1)); fprintf('\n');
end
