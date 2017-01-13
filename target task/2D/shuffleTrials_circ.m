function[T] = shuffleTrials_circ(Trials,comp_tt,trial_inds,CROSSNUM)

T.trialnums = 1:length(Trials.idx);
T.shuff_trials = randperm(length(Trials.idx));
T.CVranges=ceil(linspace(1,length(Trials.idx),CROSSNUM+1));
T.trial_inds = trial_inds;
T.lowinds = find(comp_tt(:,2)==max(comp_tt(:,2))); %% Low unc (High Kappa)
T.highinds = find(comp_tt(:,2)==min(comp_tt(:,2))); %% High unc (Low Kappa)

end