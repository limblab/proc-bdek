brain_area = 'PMd';

loop_alignments = {'target','go'};
loop_ranges = {200:200:800 , 0:200:600};

CO_index = 1;
prediction_day_indices = 2;
[percdiff, xsfp, diffb, tds, totd, tune_align, tune_ranges, neurons] = deal(cell(1,2));

for loopthrough = 1:2
    
    tune_align = loop_alignments{loopthrough};
    tune_ranges = loop_ranges{loopthrough};
    tds{loopthrough} = diff(tune_ranges);
    
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
                [cents,tune,tune_low,tune_high] = co_tuning(co_units,[nan; alldays(CO_index).tt(1:end-1,2)],...
                alldays(CO_index).tt,t1,t2,tune_align,i);

    %             [cents,tune,tune_low,tune_high] = co_tuning(co_units,alldays(CO_index).tt(:,2),...
    %             alldays(CO_index).tt,0,200,tune_align,i);
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

            neurons{loopthrough}{ranges}{i}.tuning = wrapped_tune;
            neurons{loopthrough}{ranges}{i}.tuning_low = wrapped_low;
            neurons{loopthrough}{ranges}{i}.tuning_high = wrapped_high;
        end
    end
end
%%
[expect_counts,firing_perc,trial_counts,...
    midbound,firing_diffs,good_neurso,good_neurs,unit_counts] = deal(cell(2,1));

for loopthrough = 1:2
    tune_align = loop_alignments{loopthrough};
    tune_ranges = loop_ranges{loopthrough};
    tds{loopthrough} = diff(tune_ranges);
    time_bins = tune_ranges;
    time_align = tune_align;
    
    prediction_day = prediction_day_indices;
    day_units = eval(sprintf('alldays(%d).%s_units',prediction_day,brain_area));
    if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end

    day_pred = prediction_day;
    tt_pred = alldays(day_pred).tt;

    % find rasters for each unit
    unit_pred = cell(length(day_units),1);
    for q = 1:length(day_units)

        clc;
        fprintf('Day: %d/%d\nUnit: %d/%d\n',loopthrough,2,q,length(day_units));

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
        unit_counts{loopthrough}{q} = rast./repmat(diff(time_bins)/1000,size(rast,1),1); % in FIRING RATE
    end
end
%%
for loopthrough = 1:2
    for bin = 1:length(time_bins) - 1

        tune_range = tune_ranges(bin:bin+1);
        t1 = tune_range(1);
        t2 = tune_range(2);

        clc; fprintf('bin %d/%d\n',bin,length(time_bins)-1);
        unit_cell = struct2cell(vertcat(neurons{loopthrough}{bin}{:}));    
        tuning_array = horzcat(unit_cell{1,:})'; % in FIRING RATE

        modulations = (max(tuning_array,[],2) - min(tuning_array,[],2))'; % in FR
        good_neurso{loopthrough}{bin} = find(modulations > -1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %good_neurs{bin} = good_neurs{bin}(ismember(good_neurs{bin},oklow));
        good_neurs = good_neurso;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        colors_of_plot = {'b','g','r','c'};
        %colors_of_plot = {'c','b','m','r'};

        day_pred = prediction_day_indices(1);
        tt_pred = alldays(day_pred).tt;  

        moduls = repmat(modulations,size(tt_pred,1),1);

        day_units = eval(sprintf('alldays(1).%s_units',brain_area));
        if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end

        for i = 1:size(tt_pred,1)   
            for q = 1:length(day_units)       
                trial_counts{loopthrough}{bin}(i,q) = unit_counts{loopthrough}{q}(i,bin); % in FIRING RATE
                move_aligned_thetas = abs(circ_dist(wrapped_cents,tt_pred(i,10)));
                move_theta_index = find(move_aligned_thetas==min(move_aligned_thetas),1,'first');
                expect_counts{loopthrough}{bin}(i,q) = tuning_array(q,move_theta_index); % in FIRING RATE
            end
        end

        firing_diffs{loopthrough}{bin} = (trial_counts{loopthrough}{bin}(:,good_neurs{loopthrough}{bin}) - ...
            expect_counts{loopthrough}{bin}(:,good_neurs{loopthrough}{bin}))./moduls(:,good_neurs{loopthrough}{bin});
%         firing_perc{z}{bin} = (trial_counts{z}{bin}(:,good_neurs{bin}) - ...
%             repmat(min(tuning_array(good_neurs{bin},:),[],2)',length(tt_pred),1))./moduls{z}(:,good_neurs{bin});
        firing_diffs{loopthrough}{bin}(isnan(firing_diffs{loopthrough}{bin})) = 0;
        firing_diffs{loopthrough}{bin}(isinf(firing_diffs{loopthrough}{bin})) = 0;

        %FR_change{z} = mean(firing_perc{z}{bin},2);

        set_of_inds = unique(alldays(prediction_day_indices(1)).tt(:,3));
        set_of_inds = flipud(set_of_inds);
        like_ind = cell(length(set_of_inds),1);
        for qq = 1:length(set_of_inds)

            if t1 >= 0 && t2 >= 0
                %like_ind{z}{qq} = find(tt_pred(:,3)==set_of_inds(qq));
                like_ind{qq} = find(ismember(tt_pred(:,3),set_of_inds(qq)));
            else
                like_ind{qq} = find(ismember(tt_pred(:,3),set_of_inds));
            end

            midbound{loopthrough}{qq}(:,bin) = mean(firing_diffs{loopthrough}{bin}(like_ind{qq},:),1)';  

            %totd{loopthrough}(qq,bin) = midbound{loopthrough}{qq}(:,bin);
        end
    end
end
%%
%all_mets = horzcat(midbound{:});
all_mets = horzcat(midbound{:});

cat_met = cell(size(all_mets{1},2),1);
for i = 1:length(all_mets)  
    for j = 1:size(all_mets{i},2)
        
        cat_met{j}(:,i) = all_mets{i}(:,j);
    end
end

dist_met = zeros(size(cat_met{1}));

for i = 1:(length(cat_met)-1)
    
    metric_ind = cat_met{i+1} - cat_met{1};
    dist_met = dist_met + metric_ind;
end
    
totdist_met = sum(dist_met,2);

sorted_distmet = flipud(sortrows(totdist_met));

remaining_met = zeros(length(sorted_distmet),1);
for i = 1:length(sorted_distmet)
    remaining_met(i) = sum(sorted_distmet(i:end));
end

figure; hold on;
plot((1:length(totdist_met))./length(totdist_met),remaining_met,'b.');

% plot((1:length(totdist_met))./length(totdist_met),flipud(sortrows(totdist_met)),'g.'); 
% hold on; plot([0 1],[0 0],'r--');
    
    
    





