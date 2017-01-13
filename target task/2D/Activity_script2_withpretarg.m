tune_ranges = [-800:200:800];
tune_align = 'target';
brain_area = 'M1';

neurons = cell(length(tune_ranges)-1,1);

CO_index = 1;
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

        if t1 >= 0 && t2 >= 0
            [cents,tune,tune_low,tune_high] = co_tuning(co_units,alldays(CO_index).tt(:,10),...
            alldays(CO_index).tt,t1,t2,tune_align,i);
        else
%             [cents,tune,tune_low,tune_high] = co_tuning(co_units,[nan; alldays(CO_index).tt(1:end-1,2)],...
%             alldays(CO_index).tt,t1,t2,tune_align,i);
        
            [cents,tune,tune_low,tune_high] = co_tuning(co_units,alldays(CO_index).tt(:,2),...
            alldays(CO_index).tt,0,200,tune_align,i);
        end
            
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
expect_counts = cell(length(prediction_day_indices),1);
firing_perc = cell(length(prediction_day_indices),1);
trial_counts = cell(length(prediction_day_indices),1);
lowerbound = cell(length(prediction_day_indices),1);
upperbound = cell(length(prediction_day_indices),1);
midbound = cell(length(prediction_day_indices),1);
firing_diffs = cell(length(prediction_day_indices),1);
good_neurs = cell(length(prediction_day_indices),1);

for bin = 1:length(time_bins) - 1
    clc; fprintf('bin %d/%d\n',bin,length(time_bins)-1);
    unit_cell = struct2cell(vertcat(neurons{bin}{:}));    
    tuning_array = horzcat(unit_cell{1,:})';

    modulations = (max(tuning_array,[],2) - min(tuning_array,[],2))';
    good_neurs{bin} = find(modulations > 0);

    colors_of_plot = {'b','g','r','c'};
    %figure; hold on;
    moduls = cell(length(prediction_day_indices),1);
    FR_change = moduls;
    like_ind = moduls;
    for z = 1:length(prediction_day_indices)
        
        day_pred = prediction_day_indices(z);
        tt_pred = alldays(day_pred).tt;  
        ntrials = size(tt_pred,1);
        
        tlocs = repmat(tt_pred(:,2),1,ntrials);
        averagermat = triu(ones(ntrials,ntrials),0);
        averagermat(averagermat==0) = NaN;
        tbt_prior = circ_mean(tlocs.*averagermat,[],1)';
        
        moduls{z} = repmat(modulations,size(tt_pred,1),1);

        day_units = eval(sprintf('alldays(1).%s_units',brain_area));
        if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end
 
        for i = 1:size(tt_pred,1)   
            for q = 1:length(day_units)       

                if time_bins(bin) > 0 && time_bins(bin+1) > 0
                    trial_counts{z}{bin}(i,q) = unit_counts{z}{q}(i,bin);
                    move_aligned_thetas = abs(circ_dist(wrapped_cents,tt_pred(i,10)));
                    move_theta_index = find(move_aligned_thetas==min(move_aligned_thetas),1,'first');
                    expect_counts{z}{bin}(i,q) = tuning_array(q,move_theta_index);
                else
                    trial_counts{z}{bin}(i,q) = unit_counts{z}{q}(i,bin);
                    move_aligned_thetas = abs(circ_dist(wrapped_cents,tbt_prior(i)));
                    move_theta_index = find(move_aligned_thetas==min(move_aligned_thetas),1,'first');
                    expect_counts{z}{bin}(i,q) = tuning_array(q,move_theta_index);
                end
            end
        end

        firing_diffs{z}{bin} = (trial_counts{z}{bin}(:,good_neurs{bin}) - ...
            expect_counts{z}{bin}(:,good_neurs{bin}))./moduls{z}(:,good_neurs{bin});
        firing_perc{z}{bin} = (trial_counts{z}{bin}(:,good_neurs{bin}) - ...
            repmat(min(tuning_array(good_neurs{bin},:),[],2)',length(tt_pred),1))./moduls{z}(:,good_neurs{bin});
        
        %FR_change{z} = mean(firing_diffs{z}{bin},2);
        FR_change{z} = mean(firing_perc{z}{bin},2);
    %     [kurt_low{z},kurt_high{z}] = boot_bounds(5000,@skewness,firing_diffs{z}',2.5,97.5);

        set_of_inds = unique(alldays(prediction_day_indices(z)).tt(:,3));
        set_of_inds = flipud(set_of_inds);
        for qq = 1:length(set_of_inds)
            like_ind{z}{qq} = find(tt_pred(:,3)==set_of_inds(qq));
            %[lk{z}{qq},uk{z}{qq}] = boot_bounds(10000,@mean,FR_change{z}(like_ind{z}{qq}),2.5,97.5);   
            [lowerbound{z}(qq,bin),upperbound{z}(qq,bin)] = ...
                boot_bounds(10000,@mean,FR_change{z}(like_ind{z}{qq}),2.5,97.5);     
            midbound{z}(qq,bin) = mean(FR_change{z}(like_ind{z}{qq}));
            
            %mk{z}{qq} = mean(FR_change{z}(like_ind{z}{qq}));
        end
    end
end

%%
dayind = 1;
figure; hold on;
for i = 1:length(unique(alldays(dayind+1).tt(:,3)))
    
    plot(xsforplot,midbound{dayind}(i,:),colors_of_plot{i});
    patch([xsforplot fliplr(xsforplot)],[lowerbound{dayind}(i,:) fliplr(upperbound{dayind}(i,:))],...
        colors_of_plot{i},'FaceAlpha',0.5,'EdgeAlpha',0.5);
end 
xlabel('Time From Target On (ms)','FontSize',14);
ylabel('Change in FR (% modulation)','FontSize',14);