%%
time_bins = tune_ranges;
time_align = tune_align;
prediction_day_indices = 2;
num_blocks = 7;

day_pred = prediction_day_indices;

speeds = kin_exam(alldays(day_pred).bdfM,alldays(day_pred).tt);
kinxs = 1:500;

figure; hold on;

blocks = round(linspace(1,size(alldays(prediction_day_indices).tt,1),num_blocks + 1));
%blocks = [1 25];
[midbound,lowerbound,upperbound] = deal(cell(length(blocks)-1,1));
for bs = 1:length(blocks)-1
    clc; fprintf('Block %d/%d\n',bs,length(blocks)-1);
    
    trial_block = blocks(bs):blocks(bs+1);
    tt_pred = alldays(prediction_day_indices).tt(trial_block,:); 

    day_units = eval(sprintf('alldays(1).%s_units',brain_area));
    if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end

    colors_of_plot = {'b','g','r','k'};
    
    [t_counts,e_counts,firing_diffs] = deal(cell(length(time_bins)-1,1));

    % Get Likelihood indices
    set_of_inds = unique(tt_pred(:,3));
    set_of_inds = flipud(set_of_inds);
    like_ind = cell(length(set_of_inds),1);
    for qq = 1:length(set_of_inds)
        like_ind{qq} = find(tt_pred(:,3)==set_of_inds(qq));
    end

    [lowerbound{bs},upperbound{bs},midbound{bs}] = deal(zeros(length(set_of_inds),length(time_bins)-1));
    for bin = 1:length(time_bins) - 1
        clc; fprintf('Block %d/%d: bin %d/%d\n',bs,length(blocks)-1,bin,length(time_bins)-1);
        unit_cell = struct2cell(vertcat(neurons{bin}{:}));    

        tuning_array = horzcat(unit_cell{1,:})';
        modulations = (max(tuning_array,[],2) - min(tuning_array,[],2))';
        good_neurs = find(modulations > 0);
        modls = repmat(modulations,size(tt_pred,1),1);

        % Set up trial and expected counts
        t_counts{bin} = zeros(size(tt_pred,1),length(day_units));
        e_counts{bin} = zeros(size(tt_pred,1),length(day_units));
        for i = 1:size(tt_pred,1)   
            for q = 1:length(day_units)       
                t_counts{bin}(i,q) = unit_counts{1}{q}(trial_block(i),bin);
                move_aligned_thetas = abs(circ_dist(wrapped_cents,tt_pred(i,10)));
                move_theta_index = find(move_aligned_thetas==min(move_aligned_thetas),1,'first');
                e_counts{bin}(i,q) = tuning_array(q,move_theta_index);
            end
        end

        % Differences in firing rates
        firing_diffs{bin} = (t_counts{bin}(:,good_neurs) - e_counts{bin}(:,good_neurs))./modls(:,good_neurs);        
        FR_change = mean(firing_diffs{bin},2);

        % Get boot boundaries
        for qq = 1:length(set_of_inds)
            [lowerbound{bs}(qq,bin),upperbound{bs}(qq,bin)] = boot_bounds(10000,@mean,FR_change(like_ind{qq}),2.5,97.5);     
            midbound{bs}(qq,bin) = mean(FR_change(like_ind{qq}));
        end
    end

    % Do plotting
    subplot(2,length(blocks)-1,bs); hold on;
    for i = 1:length(set_of_inds)
        plot(xsforplot,midbound{bs}(i,:),colors_of_plot{i});
        patch([xsforplot fliplr(xsforplot)],[lowerbound{bs}(i,:) fliplr(upperbound{bs}(i,:))],colors_of_plot{i},'FaceAlpha',0.5,'EdgeAlpha',0.5);
    end 
    xlabel('Time From Target On (ms)','FontSize',14);
    ylabel('Change in FR (% modulation)','FontSize',14);
    title(sprintf('%d:%d',trial_block(1),trial_block(end)),'FontSize',14);
    ylim([-0.6 1.2]);
    
    subplot(2,length(blocks)-1,bs+length(blocks)-1); hold on;
    for i = 1:length(set_of_inds)
        
        clc; fprintf('Block %d/%d: bin %d/%d: likelihood %d\n',bs,length(blocks)-1,bin,length(time_bins)-1,i);
        plot(kinxs,nanmean(speeds(like_ind{i},kinxs)),colors_of_plot{i});
        [kinlow,kinhigh] = boot_bounds(1000,@nanmean,speeds(like_ind{i},kinxs),2.5,97.5);
        
        patch([kinxs fliplr(kinxs)],[kinlow' fliplr(kinhigh')],colors_of_plot{i},'FaceAlpha',.5,'EdgeAlpha',.5);
        xlabel('Time From Go (ms)','FontSize',14);
        ylabel('Speed (cm/s)','FontSize',14);
    end

end
