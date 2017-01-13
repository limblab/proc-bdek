function[normal_tt, dropped_tt, good_indices] = dropCatchTrials(comp_tt,cutoff)

    hold_times = comp_tt(:,6) - comp_tt(:,5);
    
    good_indices = hold_times > cutoff;
    
    normal_tt = comp_tt(hold_times > cutoff,:);
    dropped_tt = comp_tt(hold_times < cutoff,:);
   
end