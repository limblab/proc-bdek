function [fit_error,model_map] = Bayes_fitting_func(Bdir,Rdir,prior_f,K_l)

d_theta = 0.01;
reps = 200;

dtc = 1/d_theta;
ts = -pi:d_theta:pi;
% Prior distribution
p_dist = prior_f(ts);
pdmat = repmat(p_dist,reps,1);

% Likelihood distribution
l_dist = circ_vmpdf(ts,0,K_l);
ldmat = zeros(length(ts));
for i = 1:length(ts)
    ldmat(i,:) = circshift(l_dist,i-round(length(ts)/2))';
end

% Preconstruct von mises lookup table
vm_lookup = zeros(1000,1);
for i = 1:1000
    vm_lookup(i) = round(dtc*vonmisrand(K_l))/dtc;
end
noise_samples = datasample(vm_lookup,reps);

theta2ind = @(x) mod(floor(dtc*x + length(ts)/2),length(ts)-1)+1;
indexing_mat = repmat(1:length(ts),reps,1);
MAP = zeros(reps,length(Bdir));
for i = 1:length(Bdir) % Loop through trials

    %theta_hats = repmat(round(dtc*Bdir(i))/dtc,reps,1) + datasample(vm_lookup,reps); % noisy estimated samples of theta
    theta_hats = repmat(ceil(dtc*Bdir(i))/dtc,reps,1) + noise_samples; 
    
    theta_hat_inds = theta2ind(theta_hats);

    posteriors = ldmat(theta_hat_inds,:).*pdmat;    
    MAP(:,i) = ts(indexing_mat(posteriors == repmat(max(posteriors,[],2),1,size(posteriors,2))))';
      
end
model_map = circ_mean(MAP,[],1)';


% fit_error = sum(circ_dist(Rdir,model_map).^2);
fit_error = circ_median(abs(circ_dist(Rdir,model_map)));
