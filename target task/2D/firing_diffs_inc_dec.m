brain_area = 'M1';

loop_alignments = {'target','go'};
loop_ranges = {[-400:200:800],[0:200:600]};

CO_index = 1;
figure; hold on;
[percdiff, xsfp, diffb] = deal(cell(1,2));
ns = zeros(2,2);
oks = cell(2,2);
for loopthrough = 1:2
    
    tune_align = loop_alignments{loopthrough};
    tune_ranges = loop_ranges{loopthrough};
    neurons = cell(length(tune_ranges)-1,1);
    
for ranges = 1:length(tune_ranges)-1
    tune_range = tune_ranges(ranges:ranges+1);
    % Initialize
    co_units = eval(sprintf('alldays(%d).%s_units',CO_index,brain_area));
    if strcmp(brain_area,'M1'); co_units(end-1:end) =[]; end
    %neurons = cell(length(co_units),length(tune_range)-1);

    for i = 1:length(co_units) % Loop through units
         clc; fprintf('%d/%d Tuning: (%d/%d)\n',ranges,length(tune_ranges)-1,i,length(co_units));

        % set time bin edges
        t1 = tune_range(1);
        t2 = tune_range(2);

        [cents,tune,tune_low,tune_high] = co_tuning(co_units,alldays(CO_index).tt(:,10),...
            alldays(CO_index).tt,t1,t2,tune_align,i);

        % Smooth and fix edges on tuning curves
        front_cent = cents;
        back_cent = cents;

        front_tune = tune; front_low = tune_low; front_high = tune_high;
        back_tune = tune; back_low = tune_low; back_high = tune_high;

        tune_cents = [back_cent-2*pi cents front_cent+2*pi];
        tune_pad = [back_tune;tune;front_tune];
        low_pad = [back_low;tune_low;front_low]; 
        high_pad = [back_high;tune_high;front_high];

        interp_cents = tune_cents(1):.01:tune_cents(end);
        tune_interp = interp1(tune_cents,tune_pad,interp_cents);
        low_interp = interp1(tune_cents,low_pad,interp_cents);
        high_interp = interp1(tune_cents,high_pad,interp_cents);

        smooth_tune = smooth(tune_interp,100);
        smooth_low = smooth(low_interp,100);
        smooth_high= smooth(high_interp,100);

        wrapped_cents = interp_cents(interp_cents >=0 & interp_cents< 2*pi);
        wrapped_tune = smooth_tune(interp_cents >= 0 & interp_cents <2*pi); 
        wrapped_low = smooth_low(interp_cents >= 0 & interp_cents <2*pi); 
        wrapped_high = smooth_high(interp_cents >= 0 & interp_cents <2*pi);

        neurons{ranges}{i}.tuning = wrapped_tune;
        neurons{ranges}{i}.tuning_low = wrapped_low;
        neurons{ranges}{i}.tuning_high = wrapped_high;
    end
end
centers = round(wrapped_cents*1000)./1000;
%%
time_bins = tune_ranges;
time_align = tune_align;
prediction_day_indices = [2];
 
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
expect_counts = cell(length(prediction_day_indices),1);
trial_counts = cell(length(prediction_day_indices),1);
[lowerbound, lowerbound1, lowerbound2] = deal(cell(length(prediction_day_indices),1));
[upperbound, upperbound1, upperbound2] = deal(cell(length(prediction_day_indices),1));
[midbound, midbound1, midbound2] = deal(cell(length(prediction_day_indices),1));
[firing_diffs, firing_diffs1, firing_diffs2] = deal(cell(length(prediction_day_indices),1));
[good_neurso, good_neurs1, good_neurs2] = deal(cell(length(prediction_day_indices),1));

for bin = 1:length(time_bins) - 1
    clc; fprintf('bin %d/%d\n',bin,length(time_bins)-1);
    unit_cell = struct2cell(vertcat(neurons{bin}{:}));    
    tuning_array = horzcat(unit_cell{1,:})';

    modulations = (max(tuning_array,[],2) - min(tuning_array,[],2))';
    good_neurso{bin} = find(modulations > 0);

    colors_of_plot = {'b','g','r','c'};
    %colors_of_plot = {'c','b','m','r'};
    moduls = cell(length(prediction_day_indices),1);
    like_ind = moduls;
    for z = 1:length(prediction_day_indices)
       
        day_pred = prediction_day_indices(z);
        tt_pred = alldays(day_pred).tt;  
       
        moduls{z} = repmat(modulations,size(tt_pred,1),1);

        day_units = eval(sprintf('alldays(1).%s_units',brain_area));
        if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end
 
        for i = 1:size(tt_pred,1)   
            for q = 1:length(day_units)       

                trial_counts{z}{bin}(i,q) = unit_counts{z}{q}(i,bin);
                move_aligned_thetas = abs(circ_dist(wrapped_cents,tt_pred(i,10)));
                move_theta_index = find(move_aligned_thetas==min(move_aligned_thetas),1,'first');
                expect_counts{z}{bin}(i,q) = tuning_array(q,move_theta_index);

            end
        end

        firing_diffs{z}{bin} = (trial_counts{z}{bin}(:,good_neurso{bin}) - ...
            expect_counts{z}{bin}(:,good_neurso{bin}))./moduls{z}(:,good_neurso{bin});
    end
end
xsforplot = .5*(time_bins(1:end-1)+time_bins(2:end));
significant_diffs_in_FR;

[FR_change, FR_change1, FR_change2] = deal(cell(length(prediction_day_indices),1));
for bin = 1:length(time_bins) - 1
    
    clc; fprintf('bin %d/%d\n',bin,length(time_bins)-1);
    unit_cell = struct2cell(vertcat(neurons{bin}{:}));    
    tuning_array = horzcat(unit_cell{1,:})';

    modulations = (max(tuning_array,[],2) - min(tuning_array,[],2))';
    moduls = cell(length(prediction_day_indices),1);
    like_ind = moduls;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    good_neurs1{bin} = good_neurso{bin}(ismember(good_neurso{bin},okhigh));
    good_neurs2{bin} = good_neurso{bin}(ismember(good_neurso{bin},oklow));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ns(:,loopthrough) = [length(okhigh); length(oklow)];
    oks{1,loopthrough} = okhigh; oks{2,loopthrough} = oklow;
    
    for z = 1:length(prediction_day_indices)
        
        day_pred = prediction_day_indices(z);
        tt_pred = alldays(day_pred).tt;  
       
        moduls{z} = repmat(modulations,size(tt_pred,1),1);

        day_units = eval(sprintf('alldays(1).%s_units',brain_area));
        if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end

        firing_diffs1{z}{bin} = (trial_counts{z}{bin}(:,good_neurs1{bin}) - ...
            expect_counts{z}{bin}(:,good_neurs1{bin}))./moduls{z}(:,good_neurs1{bin});
        
        firing_diffs2{z}{bin} = (trial_counts{z}{bin}(:,good_neurs2{bin}) - ...
            expect_counts{z}{bin}(:,good_neurs2{bin}))./moduls{z}(:,good_neurs2{bin});
        
        FR_change{z} = mean(firing_diffs{z}{bin},2);
        FR_change1{z} = mean(firing_diffs1{z}{bin},2);
        FR_change2{z} = mean(firing_diffs2{z}{bin},2);

        set_of_inds = unique(alldays(prediction_day_indices(z)).tt(:,3));
        set_of_inds = flipud(set_of_inds);
        for qq = 1:length(set_of_inds)
            like_ind{z}{qq} = find(tt_pred(:,3)==set_of_inds(qq));
            %[lk{z}{qq},uk{z}{qq}] = boot_bounds(10000,@mean,FR_change{z}(like_ind{z}{qq}),2.5,97.5);   
            [lowerbound{z}(qq,bin),upperbound{z}(qq,bin)] = ...
                boot_bounds(1000,@mean,FR_change{z}(like_ind{z}{qq}),2.5,97.5);     
            midbound{z}(qq,bin) = mean(FR_change{z}(like_ind{z}{qq}));
            
            [lowerbound1{z}(qq,bin),upperbound1{z}(qq,bin)] = ...
                boot_bounds(1000,@mean,FR_change1{z}(like_ind{z}{qq}),2.5,97.5);     
            midbound1{z}(qq,bin) = mean(FR_change1{z}(like_ind{z}{qq}));
            
            [lowerbound2{z}(qq,bin),upperbound2{z}(qq,bin)] = ...
                boot_bounds(1000,@mean,FR_change2{z}(like_ind{z}{qq}),2.5,97.5);     
            midbound2{z}(qq,bin) = mean(FR_change2{z}(like_ind{z}{qq}));
            
        end
    end
end

%%
dayind = 1;

if loopthrough == 1

    xsforplotT = .5*(time_bins(1:end-1)+time_bins(2:end));
    max_targ_x = max(xsforplotT);
    for i = 1:length(unique(alldays(dayind+1).tt(:,3)))

        subplot(1,2,1); hold on;

        plot(xsforplotT,midbound1{dayind}(i,:),colors_of_plot{i});
        patch([xsforplotT fliplr(xsforplotT)],[lowerbound1{dayind}(i,:) fliplr(upperbound1{dayind}(i,:))],...
            colors_of_plot{i},'FaceAlpha',0.5,'EdgeAlpha',0.5);

        subplot(1,2,2); hold on;
        
        plot(xsforplotT,midbound2{dayind}(i,:),colors_of_plot{i});
        patch([xsforplotT fliplr(xsforplotT)],[lowerbound2{dayind}(i,:) fliplr(upperbound2{dayind}(i,:))],...
            colors_of_plot{i},'FaceAlpha',0.5,'EdgeAlpha',0.5);
    end
    xsfp{1} = xsforplotT;

else

    xsatgo = .5*(time_bins(1:end-1)+time_bins(2:end));
    xsforplotG = xsatgo + max_targ_x + 100;
    min_go_x = min(xsatgo);

    for i = 1:length(unique(alldays(dayind+1).tt(:,3)))

        subplot(1,2,1); hold on;
        
        plot(xsforplotG,midbound1{dayind}(i,:),colors_of_plot{i});
        patch([xsforplotG fliplr(xsforplotG)],[lowerbound1{dayind}(i,:) fliplr(upperbound1{dayind}(i,:))],...
            colors_of_plot{i},'FaceAlpha',0.5,'EdgeAlpha',0.5);
        
        subplot(1,2,2); hold on;
        
        plot(xsforplotG,midbound2{dayind}(i,:),colors_of_plot{i});
        patch([xsforplotG fliplr(xsforplotG)],[lowerbound2{dayind}(i,:) fliplr(upperbound2{dayind}(i,:))],...
            colors_of_plot{i},'FaceAlpha',0.5,'EdgeAlpha',0.5);
    end 
    xsfp{2} = xsforplotG;
end

end
subplot(1,2,1);
ylim([-0.6 1.5]);
plot(max_targ_x*[1 1],[-0.6 1.5],'k--');
plot([0 0],[-0.6 1.5],'k--');
plot((max_targ_x+min_go_x+100)*[1 1],[-0.5 1.5],'k--');
xlabel('Time From Target On (ms)','FontSize',14);
ylabel('Change in FR (% modulation)','FontSize',14);
title(sprintf('%s %s increase ([%d,%d]/%d)',monkey,brain_area,ns(1,1),ns(1,2),length(modulations)),'FontSize',24);

subplot(1,2,2);
ylim([-0.6 1.5]);
plot(max_targ_x*[1 1],[-0.6 1.5],'k--');
plot([0 0],[-0.6 1.5],'k--');
plot((max_targ_x+min_go_x+100)*[1 1],[-0.5 1.5],'k--');
xlabel('Time From Target On (ms)','FontSize',14);
ylabel('Change in FR (% modulation)','FontSize',14);
title(sprintf('%s %s decrease ([%d,%d]/%d)',monkey,brain_area,ns(2,1),ns(2,2),length(modulations)),'FontSize',24);

% 
% figure; hold on;
% bar(xsforplotT,percdiff{1}); bar(xsforplotG,percdiff{2});
% plot(max_targ_x*[1 1],[0 1],'k--');
% plot([0 0],[0 1],'k--');
% xlabel('bin'); ylabel('% Different (H - L)'); 

