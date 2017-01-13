[fdif,fdif_H,fdif_L,lfdif,lfdif_L,lfdif_H] = deal(zeros(length(PMd_units),length(firing_diffs{1})));

H_inds = find(alldays(2).tt(:,3) == min(alldays(2).tt(:,3)));
L_inds = find(alldays(2).tt(:,3) == max(alldays(2).tt(:,3)));

for i = 1:length(firing_diffs{1})
    
    firediffs = firing_diffs{1}{i};
    
    fdif(good_neurso{i},i) = mean(firing_diffs{1}{i}(H_inds,:),1)';
    lfdif(good_neurso{i},i) = mean(firing_diffs{1}{i}(L_inds,:),1)';
    
    clc; fprintf('bin: %d/%d\n',i,length(firing_diffs{1}));
    [fdif_L(good_neurso{i},i), fdif_H(good_neurso{i},i)] = boot_bounds(2000,@mean,firing_diffs{1}{i}(H_inds,:),2.5,97.5);
    [lfdif_L(good_neurso{i},i), lfdif_H(good_neurso{i},i)] = boot_bounds(2000,@mean,firing_diffs{1}{i}(L_inds,:),2.5,97.5);
end

nozerbound = ~((fdif_H > 0) & (fdif_L < 0));
diffbound1 = ~((fdif_H >= lfdif) & (fdif_L <= lfdif));
diffbound2 = ~((lfdif_H >= fdif) & (lfdif_L <= fdif));
diffbound  = (diffbound1 + diffbound2) == 2; 

highbound = (fdif_L > lfdif) & (lfdif_H < fdif);
lowbound  = (fdif_H < lfdif) & (lfdif_L > fdif);

posind = find(xsforplot >= 0);

bins_of_dif = sum(diffbound(:,posind),2);
bins_of_high = sum(highbound(:,posind),2);
bins_of_low = sum(lowbound(:,posind),2);

bestneurs = find(bins_of_dif > 0.5*length(posind));

okneurs = find(bins_of_dif > 0);
okhigh = find(bins_of_high > 0);
oklow = find(bins_of_low > 0);

perc_d_inbin = sum(diffbound)/length(bins_of_dif);
perc_nz_inbin = sum(nozerbound)/length(bins_of_dif);