[modrange,modrangeH,modrangeL] = deal(cell(length(neurons{1}),1));
for timebin = 1:length(neurons)   
    for neurnum = 1:length(neurons{timebin})       
        tc = neurons{timebin}{neurnum}.tuning;
        tch = neurons{timebin}{neurnum}.tuning_high;
        tcl = neurons{timebin}{neurnum}.tuning_low;
        modrange{neurnum}(:,timebin) = [max(tc); min(tc)];     
        modrangeH{neurnum}(:,timebin) = [max(tch); min(tch)];
        modrangeL{neurnum}(:,timebin) = [max(tcl); min(tcl)];
    end
end

[vis,ebu,lbu] = deal(zeros(length(modrange),1));
for neurnum = 1:length(modrange)
    
    ranges = modrange{neurnum}(1,:) - modrange{neurnum}(2,:);
    rangesMAX = modrangeH{neurnum}(1,:) - modrangeL{neurnum}(2,:);
    
    vis(neurnum) = ranges(2) > max(rangesMAX(1)); 
    ebu(neurnum) = ranges(3) > max(rangesMAX(1));
    lbu(neurnum) = sum(ranges(7:end) > max(rangesMAX(2:4)))>0;
    
end

vis_only = find(vis == 1 & vis+ebu+lbu ==1);
with_vis = find(vis == 1);
without_vis = find(vis ~= 1);
early_build = find(ebu);
late_build = find(lbu);
any_build = unique([early_build; late_build]);
vis_plus_build = any_build(ismember(any_build,with_vis));
novis_plus_build = any_build(ismember(any_build,without_vis));

novis_plus_ebuild = early_build(ismember(early_build,without_vis));
novis_plus_lbuild = late_build(ismember(late_build,without_vis));


%%
time_bins = tune_ranges;
time_align = tune_align;
prediction_day_indices = [2];

xsforplot = .5*(time_bins(1:end-1)+time_bins(2:end));
 
unit_counts = cell(length(prediction_day_indices),1);
for z = 1:length(prediction_day_indices)

    prediction_day = prediction_day_indices(z);
    day_units = eval(sprintf('alldays(%d).%s_units',prediction_day,brain_area));
    if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end

    day_pred = prediction_day;
    tt_pred = alldays(day_pred).tt;

    % find rasters for each unit
    unit_pred = cell(length(day_units),1);
    for q = 1:length(day_units)
        
        clc;
        fprintf('Day: %d/%d\nUnit: %d/%d\n',z,length(prediction_day_indices),q,length(day_units));

        [rast_out,rast_inds] = raster_plot(day_units{q},tt_pred,[time_bins(1)./1000 time_bins(end)./1000],...
            time_align,0,'none');

        out = nan(size(vertcat(rast_out{:})));
        for m = 1:length(rast_out)
            out(rast_inds{m},:) = rast_out{m};
        end
        out(isnan(out))=0;

        rast = zeros(size(out,1),length(time_bins)-1);
        binsizes = time_bins - time_bins(1); binsizes(1) = 1;
        for v = 1:length(time_bins)-1
            rast(:,v) = sum(out(:,binsizes(v):binsizes(v+1)),2);
        end

        unit_counts{z}{q} = rast./repmat(diff(time_bins)/1000,size(rast,1),1);
    end
end

%%
good_neurs = vis_plus_build;
%good_neurs = early_build;
%good_neurs = with_vis;
 
 %(randperm(length(vis_plus_build),length(novis_plus_build)));

expect_counts = cell(length(prediction_day_indices),1);
firing_perc = cell(length(prediction_day_indices),1);
trial_counts = cell(length(prediction_day_indices),1);
lowerbound = cell(length(prediction_day_indices),1);
upperbound = cell(length(prediction_day_indices),1);
midbound = cell(length(prediction_day_indices),1);
firing_diffs = cell(length(prediction_day_indices),1);

for bin = 1:length(time_bins) - 1
    clc; fprintf('bin %d/%d\n',bin,length(time_bins)-1);
    unit_cell = struct2cell(vertcat(neurons{bin}{:}));    
    tuning_array = horzcat(unit_cell{1,:})';

    modulations = (max(tuning_array,[],2) - min(tuning_array,[],2))';
    good_neurs = good_neurs(modulations(good_neurs)>0);

    colors_of_plot = {'b','g','r','c'};
    moduls = cell(length(prediction_day_indices),1);
    FR_change = moduls;
    like_ind = moduls;
    for z = 1:length(prediction_day_indices)
        
        day_pred = prediction_day_indices(z);
        tt_pred = alldays(day_pred).tt;  
        slice_pred = alldays(day_pred).slices;
       
        moduls{z} = repmat(modulations,size(tt_pred,1),1);

        day_units = eval(sprintf('alldays(1).%s_units',brain_area));
        if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end
 
        for i = 1:size(tt_pred,1)   
            for q = 1:length(day_units)       

                trial_counts{z}{bin}(i,q) = unit_counts{z}{q}(i,bin);
                move_aligned_thetas = abs(circ_dist(wrapped_cents,tt_pred(i,9)));
                move_theta_index = find(move_aligned_thetas==min(move_aligned_thetas),1,'first');
                expect_counts{z}{bin}(i,q) = tuning_array(q,move_theta_index);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 trial_counts{z}{bin}(i,q) = unit_counts{z}{q}(i,bin);
%                 high_counts = NaN(1,5);
%                 for slicex = 1:5
%                     slice_aligned_thetas = abs(circ_dist(wrapped_cents,slice_pred(i,slicex)));
%                     slice_theta_index = find(slice_aligned_thetas==min(slice_aligned_thetas),1,'first');
%                     high_counts(slicex) = tuning_array(q,slice_theta_index);
%                 end
%                 expect_counts{z}{bin}(i,q) = sum(high_counts);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                

            end
        end

        firing_diffs{z}{bin} = (trial_counts{z}{bin}(:,good_neurs) - expect_counts{z}{bin}(:,good_neurs))./moduls{z}(:,good_neurs);              
        firing_perc{z}{bin} = (trial_counts{z}{bin}(:,good_neurs) - repmat(min(tuning_array(good_neurs,:),[],2)',length(tt_pred),1))./moduls{z}(:,good_neurs);
        
        FR_change{z} = mean(firing_diffs{z}{bin},2);
        %FR_change{z} = mean(abs(firing_diffs{z}{bin}),2);
        
        set_of_inds = unique(alldays(prediction_day_indices(z)).tt(:,3));
        set_of_inds = flipud(set_of_inds);
        for qq = 1:length(set_of_inds)
            like_ind{z}{qq} = find(tt_pred(:,3)==set_of_inds(qq));
            [lowerbound{z}(qq,bin),upperbound{z}(qq,bin)] = boot_bounds(10000,@nanmean,FR_change{z}(like_ind{z}{qq}),2.5,97.5);     
            midbound{z}(qq,bin) = nanmean(FR_change{z}(like_ind{z}{qq}));

        end
    end
end

%%
dayind = 1;
figure; hold on;
for i = 1:length(unique(alldays(dayind+1).tt(:,3)))
    
    plot(xsforplot,midbound{dayind}(i,:),colors_of_plot{i});
    patch([xsforplot fliplr(xsforplot)],[lowerbound{dayind}(i,:) fliplr(upperbound{dayind}(i,:))],colors_of_plot{i},'FaceAlpha',0.5,'EdgeAlpha',0.5);
end 
xlabel('Time From Target On (ms)','FontSize',14);
ylabel('Change in FR (% modulation)','FontSize',14);

