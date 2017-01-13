function[sse] = behavior_fit_helper(alldays,k)
sess = 2;
thets = -pi:.01:pi;

[max_like] = deal(zeros(size(alldays(sess).tt,1),1));
for i = 2:size(alldays(sess).tt,1)

    mu_l = alldays(sess).tt(i,9);
    k_l = circ_kappa(alldays(sess).slices(i,:));
    k_l = k_l.*k; 
    if k_l>500; k_l = 500; end 
      
    like_dist = circ_vmpdf(thets,mu_l,k_l);
    
    % Optional prior weighting
    mu_p = circ_mean(alldays(sess).tt(1:(i-1),2));
    k_p = circ_kappa(alldays(sess).tt(1:(i-1),2));
    prior_dist = circ_vmpdf(thets,mu_p,k_p);
    
    post_dist = prior_dist.*like_dist;
    
    max_like(i) = thets(post_dist == max(post_dist));
end
    
resid = circ_dist(alldays(sess).tt(:,10),max_like);
sse = sum(resid.^2);

%%
