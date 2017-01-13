sess = 2;
thets = -pi:.01:pi;

mu_p = circ_mean(alldays(sess).tt(:,2));
k_p = circ_kappa(alldays(sess).tt(:,2));
prior_dist = circ_vmpdf(thets,mu_p,k_p);

like_conds = unique(alldays(sess).tt(:,3));
like_colors = {'r','b','k'};

% figure; hold on;
[max_like,move_p] = deal(zeros(size(alldays(sess).tt,1),1));
for i = 2:size(alldays(sess).tt,1)

    mu_l = alldays(sess).tt(i,9);
    k_l = circ_kappa(alldays(sess).slices(i,:));
    if k_l>500; k_l = 500; end 
    
    like_dist = circ_vmpdf(thets,mu_l,k_l);
    
    % Optional prior weighting
    mu_p = circ_mean(alldays(sess).tt(1:(i-1),2));
    k_p = circ_kappa(alldays(sess).tt(1:(i-1),2));
    prior_dist = circ_vmpdf(thets,mu_p,k_p);
    
    post_dist = prior_dist.*like_dist;
    
    max_like(i) = thets(post_dist == max(post_dist));
%     plot(mu_l,max_like(i),'.','Color',like_colors{alldays(sess).tt(i,3)==like_conds}); 
    
%     move_p(i) = post_dist(abs(circ_dist(alldays(sess).tt(i,10),thets))==...
%         min(abs(circ_dist(alldays(sess).tt(i,10),thets))));
    
end

%%
figure; hold on; 
like_conds = flipud(unique(alldays(2).tt(:,3)));
col2plot = {'b','r'};
for i = 1:length(like_conds)
    is = find(alldays(2).tt(:,3)==like_conds(i));
    
    plot(alldays(2).tt(is,9),alldays(2).tt(is,10),'.','Color',col2plot{i});
    
    
    [p,S] = polyfit(alldays(2).tt(is,9),alldays(2).tt(is,10),1);
    xs = [min(alldays(2).tt(is,9)):.01:max(alldays(2).tt(is,9))];

    [Y,DELTA] = polyconf(p,xs,S,'predopt','curve');
    
    plot(xs,Y,col2plot{i});
    plot(xs,Y+DELTA,'--','Color',col2plot{i});
    plot(xs,Y-DELTA,'--','Color',col2plot{i});
    
end
plot([0 pi],[0 pi],'k--');
axis equal;
xlim([0 pi]); 
ylim([0 pi]); 