
for i = 1:length(alldays(1).PMd_units)
    tunecurve = neurons{5}{i}.tuning;
    modulation = max(tunecurve) - min(tunecurve);
    pd(i) = wrapped_cents(find(tunecurve==max(tunecurve),1));
    dfrompd = abs(angle_diff(wrapped_cents , pd(i)));
    locind(i) = find(dfrompd == min(dfrompd),1,'first');
    xloc = wrapped_cents(locind(i));
    
    centered_curve(:,i) = ([tunecurve(locind(i):end);tunecurve(1:locind(i)-1)] - min(tunecurve))/modulation;
    
    
    for j = 1:size(alldays(2).tt,1)
        
        deltas = abs(wrapped_cents - alldays(2).tt(j,10));
        tune_ind = find(deltas == min(deltas));
        fireatmove(j,i) = tunecurve(tune_ind);
        normfr(j,i) = (fireatmove(j,i) - min(tunecurve))./modulation;
        
        exp_cnt(j,i) = (expect_counts{1}{5}(j,i) - min(tunecurve))./modulation;
        act_cnt(j,i) = (trial_counts{1}{5}(j,i) - min(tunecurve))./modulation;
        
    end
end

%%
dfs = firing_diffs{1};
hfds = fds(like_ind{1}{3},:);
mfds = fds(like_ind{1}{2},:);
lfds = fds(like_ind{1}{1},:);

hindfr = mean(normfr(like_ind{1}{3},good_neurs));
mindfr = mean(normfr(like_ind{1}{2},good_neurs));
lindfr = mean(normfr(like_ind{1}{1},good_neurs));

[~,tcsorth] = sortrows(hindfr');
[~,tcsortm] = sortrows(mindfr');
[~,tcsortl] = sortrows(lindfr');

hinddfs = mean(hfds);
minddfs = mean(mfds);
linddfs = mean(lfds);

%figure; plot(hindfr,hinddfs,'b.');

figure; hold on; 
plot(bin_array(hinddfs(tcsorth),1,6),'r');
plot(bin_array(minddfs(tcsortm),1,6),'g');
plot(bin_array(linddfs(tcsortl),1,6),'b');
