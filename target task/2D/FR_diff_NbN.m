function[fdif, fdif_L, fdif_H] = FR_diff_NbN(tt,firing_diffs,units,good_neurs)

% Find indices
set_of_inds = unique(tt(:,3));
set_of_inds = flipud(set_of_inds);

% Initialize
[fdif,fdif_L,fdif_H] = deal(cell(length(set_of_inds),1));
like_ind = cell(length(set_of_inds),1);
 
for z = 1:length(set_of_inds) % loop through likelihood conditions
    
    like_ind{z} = find(tt(:,3)==set_of_inds(z)); % Get the trial indices
    
    [fdif{z},fdif_L,fdif_H] = ...
        deal(zeros(length(units),length(firing_diffs))); %Initialize outputs
    
    for i = 1:length(firing_diffs)

        % Calculate outputs
        fdif{z}(good_neurs{i},i) = mean(firing_diffs{i}(like_ind{z},:),1)';
        [fdif_L(good_neurs{i},i), fdif_H(good_neurs{i},i)] = boot_bounds(2000,@mean,firing_diffs{i}(like_ind{z},:),2.5,97.5);
    
    end
end