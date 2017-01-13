function[Est] = Bayes_ifneeded(tt,slices,tt_prev,varargin)

get_mu3_k3 = @(mu1,mu2,k1,k2) [mu1 + atan2(sin(mu2-mu1),(k1/k2)+cos(mu2-mu1)), ...
                               sqrt(k1^2 + k2^2 + 2*k1*k2*cos(mu2-mu1))];                       
prodvm = @(k1,k2,theta,mu3_k3) besseli(0,mu3_k3(2))/(2*pi*besseli(0,k1)*...
                               besseli(0,k2))*circ_vmpdf(theta,mu3_k3(1),mu3_k3(2));
                           
% variables                     
thetas = 0:0.01:2*pi;
nt = size(tt,1);

% initialize
likes = flipud(unique(tt(:,3)));
posterior = NaN(nt,length(thetas));
[Bayes_est, post_k] = deal(NaN(nt,1));
pflag = zeros(nt,1);
mu_likelihood = zeros(size(tt,1),1);

k_prior = circ_kappa(tt(:,2)); 
mu_prior = circ_mean(tt(:,2));
for i = 1:nt

    %%%%%%%%%%%%%%%%%%% LIKELIHOOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lines = slices(i,:);
    
    k_likelihood = circ_kappa(lines);
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
