%% Get Preferred Directions and non-parameterized tuning curves for tracked neurons
% Specify time bin edges used for tuning curves

tune_range = [500 750];
tune_align = 'target';
brain_area = 'PMd';

% Initialize
co_units = eval(sprintf('alldays(1).%s_units',brain_area));
if strcmp(brain_area,'M1'); co_units(end-1:end) =[]; end
neurons = cell(length(co_units),length(tune_range)-1);

for i = 1:length(co_units) % Loop through units
     clc; fprintf('Tuning: (%d/%d)\n',i,length(co_units));
     
    % set time bin edges
    t1 = tune_range(1);
    t2 = tune_range(2);

    % Run fast_pref_dir.m to calculate tuning curves. 
%     [ignore1,cents,tune,tune_se,tune_rast] = ...
%         fast_pref_dir(co_units,alldays(1).tt(:,10),...
%         alldays(1).tt,16,t1,t2,tune_align,i);

    [cents,tune,tune_low,tune_high] = co_tuning(co_units,alldays(1).tt(:,10),alldays(1).tt,t1,t2,tune_align,i);
    
    % Smooth and fix edges on tuning curves
    front_cent = cents;
    back_cent = cents;

    front_tune = tune;
    front_low = tune_low;
    front_high = tune_high;
    
    back_tune = tune;
    back_low = tune_low;
    back_high = tune_high;
    
%     front_tune = repmat(tune(1),1,3); 
%     front_low = repmat(tune_low(1),1,3);
%     front_high = repmat(tune_high(1),1,3);
% 
%     back_tune = repmat(tune(end),1,3);
%     back_low = repmat(tune_low(end),1,3);
%     back_high = repmat(tune_high(end),1,3);

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

    neurons{i}.tuning = wrapped_tune;
    neurons{i}.tuning_low = wrapped_low;
    neurons{i}.tuning_high = wrapped_high;
    %neurons{i}.raster = tune_rast;
   
end
centers = round(wrapped_cents*1000)./1000;


%% Get Spike Counts
bins = {-800:100:1000 , 0:100:1000};
aligns = {'target','go'};
prediction_day_indices = [2,3];
prior_mean = 90;

confs = cell(length(bins),1);
ee = cell(length(bins),1);
est_dirs = cell(length(bins),1);
stdevs = cell(length(bins),1);
probs = cell(length(bins),1);
max_like = cell(length(bins),1);
pos_like = cell(length(bins),1);
pri_like = cell(length(bins),1);
pri_pref = cell(length(bins),1);
prior_mean = prior_mean/180*pi;
for jj = 1:length(bins)

    time_bins = bins{jj};
    time_align = aligns{jj};
    
    DAYS_pred = cell(length(prediction_day_indices),1);
    for z = 1:length(prediction_day_indices)

        prediction_day = prediction_day_indices(z);
        day_units = eval(sprintf('alldays(%d).%s_units',prediction_day,brain_area));
        if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end

        day_pred = prediction_day;
        tt_pred = alldays(day_pred).tt;

        % find rasters for each unit
        unit_pred = cell(length(day_units),1);
        for q = 1:length(day_units)

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

            unit_pred{q} = rast;
            DAYS_pred{z}{q} = rast;

            clc; fprintf('Spike Counts - %s:session %d (%d/%d)\n',aligns{jj},z,q,length(day_units));

        end
    end


    % Do population code estimates
    DAYS = prediction_day_indices;
    booted = cell(length(DAYS),1);
    conds = cell(length(DAYS),1);
    conds_cos = cell(length(DAYS),1);
    cond_sum = cell(length(DAYS),1);
    estimated_dir = cell(length(DAYS),1);
    actual_dir = cell(length(DAYS),1);
    actual_cen = cell(length(DAYS),1);
    estimate_error_dir = cell(length(DAYS),1);
    estimate_error_cen = cell(length(DAYS),1);
    estimate_conf = cell(length(DAYS),1);
    pred_conf = cell(length(DAYS),1);
    pred_stdevs = cell(length(DAYS),1);
    prob_of_pos = cell(length(DAYS),1);
    prob_of_pri = cell(length(DAYS),1);
    pref_of_pri = cell(length(DAYS),1);
    for daynumber = 1:length(DAYS)

        prediction_day = DAYS(daynumber);

        day_tt = alldays(prediction_day).tt;
        if length(unique(day_tt(:,3)))==1
            day_tt = repmat(day_tt,length(unique(day_tt(:,3))),1);
        end
        
        %Ps = exp(prior_kappa.*cos(prior_mean-wrapped_cents))./(2.*pi.*besseli(0,prior_kappa));
        Ps = ones(1,length(wrapped_cents));
        %Ps(wrapped_cents > pi) = 0;
        conds{daynumber} = cell(length(time_bins)-1,1);
        cond_sum{daynumber} = cell(length(time_bins)-1,1);
        conds_cos{daynumber} = cell(length(time_bins)-1,1);
        estimated_dir{daynumber} = zeros(length(day_tt),length(time_bins)-1);
        estimate_error_dir{daynumber} = zeros(length(day_tt),length(time_bins)-1);
        estimate_error_cen{daynumber} = zeros(length(day_tt),length(time_bins)-1);
        estimate_conf{daynumber} = zeros(length(day_tt),length(time_bins)-1);
        actual_dir{daynumber} = zeros(length(day_tt),1);
        actual_cen{daynumber} = zeros(length(day_tt),1);
        pred_conf{daynumber} = zeros(length(day_tt),1);
        pred_stdevs{daynumber} = zeros(length(day_tt),1);
        prob_of_pos{daynumber} = zeros(length(day_tt),1);
        prob_of_pri{daynumber} = zeros(length(day_tt),1);
        pref_of_pri{daynumber} = zeros(length(day_tt),1);

        dt = diff(time_bins(1:2))./1000;

        num_reps = 10000;
        unit_pred = DAYS_pred{daynumber};
        for bin_num = 1:length(time_bins)-1;

            clc; fprintf('Estimates: session %d (%d/%d)\n',daynumber,bin_num,length(time_bins)-1);

            unit_cell = struct2cell(vertcat(neurons{:,1}));    
            tuning_array = horzcat(unit_cell{1,:})';

            spike_counts = zeros(length(unit_pred),length(day_tt));
            for i = 1:length(unit_pred)
                spike_counts(i,:) = unit_pred{i}(:,bin_num)';
            end

            for i = 1:length(day_tt)
                N = diff(wrapped_cents(1:2));

    %             tuning_array = tuning_array(1:20,:);
    %             spike_counts = spike_counts(1:20,:);

                paren = Ps.*prod(exp(-tuning_array.*dt),1);
                ta_pow = tuning_array.^repmat(spike_counts(:,i),1,length(wrapped_cents));
                ta_pow = ta_pow./repmat(max(ta_pow,[],2),1,size(ta_pow,2));
                ta_pow(isnan(ta_pow))=1;
                calc_post = paren.*prod(ta_pow,1);
                %calc_post = Ps.*prod(exp(-tuning_array.*dt),1).*prod(tuning_array.^repmat(spike_counts(:,i),1,length(wrapped_cents)),1);

                conds{daynumber}{bin_num}(i,:) = calc_post./(sum(calc_post)*N);

                with_nans = conds{daynumber}{bin_num}(i,:);
                with_nans(isinf(with_nans))=0;
                with_nans(isnan(with_nans))=0;
                samps = sample_from_pdf(wrapped_cents,with_nans,10000);

                %estimated_dir{daynumber}(i,bin_num) = mean_angle(samps,'rads');
                estimated_dir{daynumber}(i,bin_num) = wrapped_cents(find(with_nans == max(with_nans),1,'first'));

                actual_dir{daynumber}(i) = day_tt(i,10);
                actual_cen{daynumber}(i) = day_tt(i,9);

%                 estimate_error_dir{daynumber}(i,bin_num) = ...
%                     angle_diff(estimated_dir{daynumber}(i,bin_num),actual_dir{daynumber}(i));
                
                estimate_error_dir{daynumber}(i,bin_num) = ...
                circ_dist(estimated_dir{daynumber}(i,bin_num),actual_dir{daynumber}(i));

%                 estimate_error_cen{daynumber}(i,bin_num) = ...
%                     angle_diff(estimated_dir{daynumber}(i,bin_num),actual_cen{daynumber}(i));
                
                estimate_error_cen{daynumber}(i,bin_num) = ...
                    circ_dist(estimated_dir{daynumber}(i,bin_num),actual_cen{daynumber}(i));

                estimate_conf{daynumber}(i,bin_num) = max(with_nans);
                pred_conf{daynumber}(i,bin_num) = kappa_calc(samps,0.05);
                pred_stdevs{daynumber}(i,bin_num) = circ_std(samps);
                
                dists_to_end = abs(circ_dist(wrapped_cents,actual_dir{daynumber}(i)));
                dists_to_pri = abs(circ_dist(wrapped_cents,prior_mean));
                dists_to_lik = abs(circ_dist(wrapped_cents,actual_cen{daynumber}(i)));
                
                true_pos_ind = find(dists_to_end == min(dists_to_end),1,'first');
                pri_pos_ind  = find(dists_to_pri == min(dists_to_pri),1,'first');
                lik_pos_ind  = find(dists_to_lik == min(dists_to_lik),1,'first');
                
                prob_of_pos{daynumber}(i,bin_num) = with_nans(true_pos_ind);
                prob_of_pri{daynumber}(i,bin_num) = with_nans(pri_pos_ind);
                pref_of_pri{daynumber}(i,bin_num) = with_nans(pri_pos_ind)-with_nans(true_pos_ind);
            end     
        end
    end

    ee{jj} = estimate_error_dir;
    est_dirs{jj} = estimated_dir;
    confs{jj} = pred_conf;
    stdevs{jj} = pred_stdevs;
    probs{jj} = conds;
    max_like{jj} = estimate_conf;
    pos_like{jj} = prob_of_pos;
    pri_like{jj} = prob_of_pri;
    pri_pref{jj} = pref_of_pri;
    
end
clc;
%% Plot Absolute Error for two epochs
figure; hold on;
xsforplot_target = -750:100:950;
xsforplot_go = 950 + (50:100:950);
for i = 1:length(DAYS);
    
    day_tt = alldays(DAYS(i)).tt;
    
    linds = find(day_tt(:,3)==max(day_tt(:,3)));
    hinds = find(day_tt(:,3)==min(day_tt(:,3)));
    
    % TARGET epoch      
    
    [dml.t, dul_l.t, dll_l.t] = circ_mean(abs(ee{1}{i}(linds,:)));
    targ_L_std = circ_std(abs(ee{1}{i}(linds,:)));
    
    [dmh.t, dul_h.t, dll_h.t] = circ_mean(abs(ee{1}{i}(hinds,:)));
    targ_H_std = circ_std(abs(ee{1}{i}(hinds,:)));
    
    % GO epoch
    [dml.g, dul_l.g, dll_l.g] = circ_mean(abs(ee{2}{i}(linds,:)));
    go_L_std = circ_std(abs(ee{2}{i}(linds,:)));
    [dmh.g, dul_h.g, dll_h.g] = circ_mean(abs(ee{2}{i}(hinds,:)));
    go_H_std = circ_std(abs(ee{2}{i}(hinds,:)));
    
    %PLOT Target
    subplot(1,length(DAYS),i); hold on;
    patch([xsforplot_target fliplr(xsforplot_target)],[dll_l.t fliplr(dul_l.t)],'b','FaceAlpha',0.5,'EdgeAlpha',0.5);
    patch([xsforplot_target fliplr(xsforplot_target)],[dll_h.t fliplr(dul_h.t)],'r','FaceAlpha',0.5,'EdgeAlpha',0.5);
% 
%     patch([xsforplot_target fliplr(xsforplot_target)],[dml.t+1.95*targ_L_std fliplr(dml.t-1.95*targ_L_std)],'b','FaceAlpha',0.5,'EdgeAlpha',0.5);
%     patch([xsforplot_target fliplr(xsforplot_target)],[dml.t+1.95*targ_H_std fliplr(dml.t-1.95*targ_H_std)],'r','FaceAlpha',0.5,'EdgeAlpha',0.5);

    plot(xsforplot_target,dml.t,'b');
    plot(xsforplot_target,dmh.t,'r');
    
    %PLOT Go
    patch([xsforplot_go fliplr(xsforplot_go)],[dll_l.g fliplr(dul_l.g)],'b','FaceAlpha',0.5,'EdgeAlpha',0.5);
    patch([xsforplot_go fliplr(xsforplot_go)],[dll_h.g fliplr(dul_h.g)],'r','FaceAlpha',0.5,'EdgeAlpha',0.5);
    
%     patch([xsforplot_go fliplr(xsforplot_go)],[dml.g+1.95*go_L_std fliplr(dml.g-1.95*go_L_std)],'b','FaceAlpha',0.5,'EdgeAlpha',0.5);
%     patch([xsforplot_go fliplr(xsforplot_go)],[dml.g+1.95*go_H_std fliplr(dml.g-1.95*go_H_std)],'r','FaceAlpha',0.5,'EdgeAlpha',0.5);
    
    plot(xsforplot_go,dml.g,'b');
    plot(xsforplot_go,dmh.g,'r');
    plot([950 950],[0 pi],'k--');
    plot([1000 1000],[0 pi],'k--');
    
    title(sprintf('End Position Error (Prior: %d)',priors{DAYS(i)}.val),'FontSize',16);
    xlabel('Time from Target (ms)','FontSize',14);
    ylabel('abs(prediction error)','FontSize',14);
    plot([0 0],[0 pi],'g--');
    ylim([0 pi]);

end

%% Plot Error for two epochs
figure; hold on;
xsforplot_target = -750:100:950;
xsforplot_go = 950 + (50:100:950);
for i = 1:length(DAYS);
    
    day_tt = alldays(DAYS(i)).tt;
    
    linds = find(day_tt(:,3)==max(day_tt(:,3)));
    hinds = find(day_tt(:,3)==min(day_tt(:,3)));
    
    % TARGET epoch      
    
    [dml.t, dul_l.t, dll_l.t] = circ_mean(ee{1}{i}(linds,:));
    targ_L_std = circ_std(ee{1}{i}(linds,:));
    [dmh.t, dul_h.t, dll_h.t] = circ_mean(ee{1}{i}(hinds,:));
    targ_H_std = circ_std(ee{1}{i}(hinds,:));
    
    % GO epoch
    [dml.g, dul_l.g, dll_l.g] = circ_mean(ee{2}{i}(linds,:));
    go_L_std = circ_std(ee{2}{i}(linds,:));
    [dmh.g, dul_h.g, dll_h.g] = circ_mean(ee{2}{i}(hinds,:));
    go_H_std = circ_std(ee{2}{i}(hinds,:));
    
    %PLOT Target
    subplot(1,length(DAYS),i); hold on;
    patch([xsforplot_target fliplr(xsforplot_target)],[dll_l.t fliplr(dul_l.t)],'b','FaceAlpha',0.5,'EdgeAlpha',0.5);
    patch([xsforplot_target fliplr(xsforplot_target)],[dll_h.t fliplr(dul_h.t)],'r','FaceAlpha',0.5,'EdgeAlpha',0.5);
% 
%     patch([xsforplot_target fliplr(xsforplot_target)],[dml.t+1.95*targ_L_std fliplr(dml.t-1.95*targ_L_std)],'b','FaceAlpha',0.5,'EdgeAlpha',0.5);
%     patch([xsforplot_target fliplr(xsforplot_target)],[dml.t+1.95*targ_H_std fliplr(dml.t-1.95*targ_H_std)],'r','FaceAlpha',0.5,'EdgeAlpha',0.5);

    plot(xsforplot_target,dml.t,'b');
    plot(xsforplot_target,dmh.t,'r');
    
    %PLOT Go
    patch([xsforplot_go fliplr(xsforplot_go)],[dll_l.g fliplr(dul_l.g)],'b','FaceAlpha',0.5,'EdgeAlpha',0.5);
    patch([xsforplot_go fliplr(xsforplot_go)],[dll_h.g fliplr(dul_h.g)],'r','FaceAlpha',0.5,'EdgeAlpha',0.5);
    
%     patch([xsforplot_go fliplr(xsforplot_go)],[dml.g+1.95*go_L_std fliplr(dml.g-1.95*go_L_std)],'b','FaceAlpha',0.5,'EdgeAlpha',0.5);
%     patch([xsforplot_go fliplr(xsforplot_go)],[dml.g+1.95*go_H_std fliplr(dml.g-1.95*go_H_std)],'r','FaceAlpha',0.5,'EdgeAlpha',0.5);
    
    plot(xsforplot_go,dml.g,'b');
    plot(xsforplot_go,dmh.g,'r');
    plot([950 950],[-5 5],'k--');
    plot([1000 1000],[-5 5],'k--');
    
    title(sprintf('End Position Error (Prior: %d)',priors{DAYS(i)}.val),'FontSize',16);
    xlabel('Time from Target (ms)','FontSize',14);
    ylabel('prediction error','FontSize',14);
    plot([0 0],[-5 5],'g--');
    ylim([-5 5]);

end
%% Concentration
figure;
for j = 1:2
    linds = find(alldays(prediction_day_indices(j)).tt(:,3) == max(alldays(prediction_day_indices(j)).tt(:,3)));
    hinds = find(alldays(prediction_day_indices(j)).tt(:,3) == min(alldays(prediction_day_indices(j)).tt(:,3)));

    for i = 1:size(est_dirs{1}{1},2)
        [tk_LOW(j,i), tk_LOW_l(j,i), tk_LOW_h(j,i)] = kappa_calc(ee{1}{j}(linds,i),0.05);
        [tk_HIGH(j,i), tk_HIGH_l(j,i), tk_HIGH_h(j,i)] = kappa_calc(ee{1}{j}(hinds,i),0.05);
    end

    for i = 1:size(est_dirs{2}{1},2)
        [gk_LOW(j,i), gk_LOW_l(j,i), gk_LOW_h(j,i)] = kappa_calc(ee{2}{j}(linds,i),0.05);
        [gk_HIGH(j,i), gk_HIGH_l(j,i), gk_HIGH_h(j,i)] = kappa_calc(ee{2}{j}(hinds,i),0.05);
    end

end

for j = 1:2
    
    subplot(1,2,j); hold on;
    
    plot(xsforplot_target,tk_LOW(j,:),'b'); 
    patch([xsforplot_target fliplr(xsforplot_target)],[tk_LOW_l(j,:) fliplr(tk_LOW_h(j,:))],'b','FaceAlpha',0.5,'EdgeAlpha',0.5); 

    plot(xsforplot_target,tk_HIGH(j,:),'r'); 
    patch([xsforplot_target fliplr(xsforplot_target)],[tk_HIGH_l(j,:) fliplr(tk_HIGH_h(j,:))],'r','FaceAlpha',0.5,'EdgeAlpha',0.5); 

    plot(xsforplot_go,gk_LOW(j,:),'b'); 
    patch([xsforplot_go fliplr(xsforplot_go)],[gk_LOW_l(j,:) fliplr(gk_LOW_h(j,:))],'b','FaceAlpha',0.5,'EdgeAlpha',0.5); 

    plot(xsforplot_go,gk_HIGH(j,:),'r'); 
    patch([xsforplot_go fliplr(xsforplot_go)],[gk_HIGH_l(j,:) fliplr(gk_HIGH_h(j,:))],'r','FaceAlpha',0.5,'EdgeAlpha',0.5); 

    
    title(sprintf('Error Concentration (Prior: %d)',priors{DAYS(j)}.val),'FontSize',16);
    xlabel('Time from Target (ms)','FontSize',14);
    ylabel('Error Concentratoin','FontSize',14);
    plot([0 0],[0 10],'g--');
    plot([950 950],[0 10],'k--');
    plot([1000 1000],[0 10],'k--');
    ylim([0 10]);
end

%% Maximum Likelihood
figure;
for j = 1:2
    linds = find(alldays(prediction_day_indices(j)).tt(:,3) == max(alldays(prediction_day_indices(j)).tt(:,3)));
    hinds = find(alldays(prediction_day_indices(j)).tt(:,3) == min(alldays(prediction_day_indices(j)).tt(:,3)));

    for i = 1:size(est_dirs{1}{1},2)
        tml_LOW(j,i) = nanmean(max_like{1}{j}(linds,i));
        tml_LOW_se(j,i) = 1.95*nanstd(max_like{1}{j}(linds,i))./sqrt(length(linds));
        tml_HIGH(j,i)= nanmean(max_like{1}{j}(hinds,i));
        tml_HIGH_se(j,i) = 1.95*nanstd(max_like{1}{j}(hinds,i))./sqrt(length(hinds));
    end

    for i = 1:size(est_dirs{2}{1},2)
        gml_LOW(j,i) = nanmean(max_like{2}{j}(linds,i));
        gml_LOW_se(j,i) = 1.95*nanstd(max_like{2}{j}(linds,i))./sqrt(length(linds));
        gml_HIGH(j,i)= nanmean(max_like{2}{j}(hinds,i));
        gml_HIGH_se(j,i) = 1.95*nanstd(max_like{2}{j}(hinds,i))./sqrt(length(hinds));
    end

end

for j = 1:2
    
    subplot(1,2,j); hold on;
    
    plot(xsforplot_target,tml_LOW(j,:),'b'); 
    patch([xsforplot_target fliplr(xsforplot_target)],...
        [tml_LOW(j,:)-tml_LOW_se(j,:) fliplr(tml_LOW(j,:)+tml_LOW_se(j,:))],'b','FaceAlpha',0.5,'EdgeAlpha',0.5); 

    plot(xsforplot_target,tml_HIGH(j,:),'r'); 
    patch([xsforplot_target fliplr(xsforplot_target)],...
        [tml_HIGH(j,:)-tml_HIGH_se(j,:) fliplr(tml_HIGH(j,:)+tml_HIGH_se(j,:))],'r','FaceAlpha',0.5,'EdgeAlpha',0.5); 

    plot(xsforplot_go,gml_LOW(j,:),'b'); 
    patch([xsforplot_go fliplr(xsforplot_go)],...
        [gml_LOW(j,:)-gml_LOW_se(j,:) fliplr(gml_LOW(j,:)+gml_LOW_se(j,:))],'b','FaceAlpha',0.5,'EdgeAlpha',0.5); 

    plot(xsforplot_go,gml_HIGH(j,:),'r'); 
    patch([xsforplot_go fliplr(xsforplot_go)],...
        [gml_HIGH(j,:)-gml_HIGH_se(j,:) fliplr(gml_HIGH(j,:)+gml_HIGH_se(j,:))],'r','FaceAlpha',0.5,'EdgeAlpha',0.5); 

    
    title(sprintf('Maximum Likelihood (Prior: %d)',priors{DAYS(j)}.val),'FontSize',16);
    xlabel('Time from Target (ms)','FontSize',14);
    ylabel('Maximum Likelihood','FontSize',14);
    plot([0 0],[0 5],'g--');
    plot([950 950],[0 5],'k--');
    plot([1000 1000],[0 5],'k--');
    ylim([0 5]);
end

%% Prediction Kappas
figure;
for j = 1:2
    linds = find(alldays(prediction_day_indices(j)).tt(:,3) == max(alldays(prediction_day_indices(j)).tt(:,3)));
    hinds = find(alldays(prediction_day_indices(j)).tt(:,3) == min(alldays(prediction_day_indices(j)).tt(:,3)));

    for i = 1:size(est_dirs{1}{1},2)
        tkap_LOW(j,i) = nanmean(confs{1}{j}(linds,i));
        [tkap_LOW_lower(j,i),tkap_LOW_upper(j,i)] = boot_bounds(10000,@nanmean,confs{1}{j}(linds,i),2.5,97.5);
        tkap_HIGH(j,i)= nanmean(confs{1}{j}(hinds,i));
        [tkap_HIGH_lower(j,i),tkap_HIGH_upper(j,i)] = boot_bounds(10000,@nanmean,confs{1}{j}(hinds,i),2.5,97.5);
    end

    for i = 1:size(est_dirs{2}{1},2)
        gkap_LOW(j,i) = nanmean(confs{2}{j}(linds,i));
        [gkap_LOW_lower(j,i), gkap_LOW_upper(j,i)] = boot_bounds(10000,@nanmean,confs{2}{j}(linds,i),2.5,97.5);
        gkap_HIGH(j,i)= nanmean(confs{2}{j}(hinds,i));
        [gkap_HIGH_lower(j,i), gkap_HIGH_upper(j,i)] = boot_bounds(10000,@nanmean,confs{2}{j}(hinds,i),2.5,97.5);
    end

end

for j = 1:2
    
    subplot(1,2,j); hold on;
    
    plot(xsforplot_target,tkap_LOW(j,:),'b'); 
    patch([xsforplot_target fliplr(xsforplot_target)],...
        [tkap_LOW_lower(j,:) fliplr(tkap_LOW_upper(j,:))],'b','FaceAlpha',0.5,'EdgeAlpha',0.5); 

    plot(xsforplot_target,tkap_HIGH(j,:),'r'); 
    patch([xsforplot_target fliplr(xsforplot_target)],...
        [tkap_HIGH_lower(j,:) fliplr(tkap_HIGH_upper(j,:))],'r','FaceAlpha',0.5,'EdgeAlpha',0.5); 

    plot(xsforplot_go,gkap_LOW(j,:),'b'); 
    patch([xsforplot_go fliplr(xsforplot_go)],...
        [gkap_LOW_lower(j,:) fliplr(gkap_LOW_upper(j,:))],'b','FaceAlpha',0.5,'EdgeAlpha',0.5); 

    plot(xsforplot_go,gkap_HIGH(j,:),'r'); 
    patch([xsforplot_go fliplr(xsforplot_go)],...
        [gkap_HIGH_lower(j,:) fliplr(gkap_HIGH_upper(j,:))],'r','FaceAlpha',0.5,'EdgeAlpha',0.5); 

    
    title(sprintf('Prediction Kappa (Prior: %d)',priors{DAYS(j)}.val),'FontSize',16);
    xlabel('Time from Target (ms)','FontSize',14);
    ylabel('Prediction Kappa','FontSize',14);
    plot([0 0],[0 50],'g--');
    plot([950 950],[0 50],'k--');
    plot([1000 1000],[0 50],'k--');
    ylim([0 50]);
end

%% position probability
figure;
xsforplot_target = -750:100:950;
xsforplot_go = 950 + (50:100:950);
for j = 1:2
    set_of_likes = unique(alldays(prediction_day_indices(j)).tt(:,3));
    
    if length(set_of_likes) > 1
        linds = find(alldays(prediction_day_indices(j)).tt(:,3) == max(set_of_likes));
        hinds = find(alldays(prediction_day_indices(j)).tt(:,3) == min(set_of_likes));
        %minds = find(~ismember(alldays(prediction_day_indices(j)).tt(:,3),[max(set_of_likes) min(set_of_likes)]));
    else
      linds = 1:length(alldays(prediction_day_indices(j)).tt);
      hinds = linds;
      %minds = linds;
    end
        
    for i = 1:size(est_dirs{1}{1},2)
        tpos_LOW(j,i) = nanmean(pos_like{1}{j}(linds,i));
        [tpos_LOW_lower(j,i),tpos_LOW_upper(j,i)] = boot_bounds(10000,@nanmean,pos_like{1}{j}(linds,i),2.5,97.5);
        tpos_HIGH(j,i)= nanmean(pos_like{1}{j}(hinds,i));
        [tpos_HIGH_lower(j,i),tpos_HIGH_upper(j,i)] = boot_bounds(10000,@nanmean,pos_like{1}{j}(hinds,i),2.5,97.5);
        %%% MEDIUM LIKE STUFF
%         tpos_MED(j,i)= nanmean(pos_like{1}{j}(minds,i));
%         [tpos_MED_lower(j,i),tpos_MED_upper(j,i)] = boot_bounds(10000,@nanmean,pos_like{1}{j}(minds,i),2.5,97.5);
        
    end

    for i = 1:size(est_dirs{2}{1},2)
        gpos_LOW(j,i) = nanmean(pos_like{2}{j}(linds,i));
        [gpos_LOW_lower(j,i), gpos_LOW_upper(j,i)] = boot_bounds(10000,@nanmean,pos_like{2}{j}(linds,i),2.5,97.5);
        gpos_HIGH(j,i)= nanmean(pos_like{2}{j}(hinds,i));
        [gpos_HIGH_lower(j,i), gpos_HIGH_upper(j,i)] = boot_bounds(10000,@nanmean,pos_like{2}{j}(hinds,i),2.5,97.5);
        %%% MEDIUM LIKE STUFF
%         gpos_MED(j,i)= nanmean(pos_like{2}{j}(minds,i));
%         [gpos_MED_lower(j,i), gpos_MED_upper(j,i)] = boot_bounds(10000,@nanmean,pos_like{2}{j}(minds,i),2.5,97.5);
    end

end

for j = 1:2
    
    subplot(1,2,j); hold on;
    
    plot(xsforplot_target,tpos_LOW(j,:),'b'); 
    patch([xsforplot_target fliplr(xsforplot_target)],...
        [tpos_LOW_lower(j,:) fliplr(tpos_LOW_upper(j,:))],'b','FaceAlpha',0.5,'EdgeAlpha',0.5); 

    plot(xsforplot_target,tpos_HIGH(j,:),'r'); 
    patch([xsforplot_target fliplr(xsforplot_target)],...
        [tpos_HIGH_lower(j,:) fliplr(tpos_HIGH_upper(j,:))],'r','FaceAlpha',0.5,'EdgeAlpha',0.5); 
    
    %%% MEDIUM LIKE
%     plot(xsforplot_target,tpos_MED(j,:),'g'); 
%     patch([xsforplot_target fliplr(xsforplot_target)],...
%         [tpos_MED_lower(j,:) fliplr(tpos_MED_upper(j,:))],'g','FaceAlpha',0.5,'EdgeAlpha',0.5); 
    %%%
    
    plot(xsforplot_go,gpos_LOW(j,:),'b'); 
    patch([xsforplot_go fliplr(xsforplot_go)],...
        [gpos_LOW_lower(j,:) fliplr(gpos_LOW_upper(j,:))],'b','FaceAlpha',0.5,'EdgeAlpha',0.5); 

    plot(xsforplot_go,gpos_HIGH(j,:),'r'); 
    patch([xsforplot_go fliplr(xsforplot_go)],...
        [gpos_HIGH_lower(j,:) fliplr(gpos_HIGH_upper(j,:))],'r','FaceAlpha',0.5,'EdgeAlpha',0.5); 

    %%% MEDIUM LIKE
%     plot(xsforplot_go,gpos_MED(j,:),'g'); 
%     patch([xsforplot_go fliplr(xsforplot_go)],...
%         [gpos_MED_lower(j,:) fliplr(gpos_MED_upper(j,:))],'g','FaceAlpha',0.5,'EdgeAlpha',0.5); 
    %%%
    
    title(sprintf('Probability of Endpoint (Prior: %d)',priors{DAYS(j)}.val),'FontSize',16);
    xlabel('Time from Target (ms)','FontSize',14);
    ylabel('Probability','FontSize',14);
    plot([0 0],[0 1.2],'g--');
    plot([950 950],[0 1.2],'k--');
    plot([1000 1000],[0 1.2],'k--');
    ylim([0 1.2]);
end
%%
for i = 1:length(pos_like)
    
    for j = 1:length(pos_like{1})
        
        max_prob{i}{j} = max(pos_like{i}{j},[],2);

    end
end

%% prior probability
figure;
xsforplot_target = -750:100:950;
xsforplot_go = 950 + (50:100:950);
for j = 1:2
    set_of_likes = unique(alldays(prediction_day_indices(j)).tt(:,3));
    
    if length(set_of_likes) > 1
        linds = find(alldays(prediction_day_indices(j)).tt(:,3) == max(set_of_likes));
        hinds = find(alldays(prediction_day_indices(j)).tt(:,3) == min(set_of_likes));
        %minds = find(~ismember(alldays(prediction_day_indices(j)).tt(:,3),[max(set_of_likes) min(set_of_likes)]));
    else
      linds = 1:length(alldays(prediction_day_indices(j)).tt);
      hinds = linds;
      %minds = linds;
    end
        
    for i = 1:size(est_dirs{1}{1},2)
        tpos_LOW(j,i) = nanmean(pri_like{1}{j}(linds,i));
        [tpos_LOW_lower(j,i),tpos_LOW_upper(j,i)] = boot_bounds(10000,@nanmean,pri_like{1}{j}(linds,i),2.5,97.5);
        tpos_HIGH(j,i)= nanmean(pri_like{1}{j}(hinds,i));
        [tpos_HIGH_lower(j,i),tpos_HIGH_upper(j,i)] = boot_bounds(10000,@nanmean,pri_like{1}{j}(hinds,i),2.5,97.5);
        %%% MEDIUM LIKE STUFF
%         tpos_MED(j,i)= nanmean(pos_like{1}{j}(minds,i));
%         [tpos_MED_lower(j,i),tpos_MED_upper(j,i)] = boot_bounds(10000,@nanmean,pos_like{1}{j}(minds,i),2.5,97.5);
        
    end

    for i = 1:size(est_dirs{2}{1},2)
        gpos_LOW(j,i) = nanmean(pri_like{2}{j}(linds,i));
        [gpos_LOW_lower(j,i), gpos_LOW_upper(j,i)] = boot_bounds(10000,@nanmean,pri_like{2}{j}(linds,i),2.5,97.5);
        gpos_HIGH(j,i)= nanmean(pri_like{2}{j}(hinds,i));
        [gpos_HIGH_lower(j,i), gpos_HIGH_upper(j,i)] = boot_bounds(10000,@nanmean,pri_like{2}{j}(hinds,i),2.5,97.5);
        %%% MEDIUM LIKE STUFF
%         gpos_MED(j,i)= nanmean(pos_like{2}{j}(minds,i));
%         [gpos_MED_lower(j,i), gpos_MED_upper(j,i)] = boot_bounds(10000,@nanmean,pos_like{2}{j}(minds,i),2.5,97.5);
    end

end

for j = 1:2
    
    subplot(1,2,j); hold on;
    
    plot(xsforplot_target,tpos_LOW(j,:),'b'); 
    patch([xsforplot_target fliplr(xsforplot_target)],...
        [tpos_LOW_lower(j,:) fliplr(tpos_LOW_upper(j,:))],'b','FaceAlpha',0.5,'EdgeAlpha',0.5); 

    plot(xsforplot_target,tpos_HIGH(j,:),'r'); 
    patch([xsforplot_target fliplr(xsforplot_target)],...
        [tpos_HIGH_lower(j,:) fliplr(tpos_HIGH_upper(j,:))],'r','FaceAlpha',0.5,'EdgeAlpha',0.5); 
    
    %%% MEDIUM LIKE
%     plot(xsforplot_target,tpos_MED(j,:),'g'); 
%     patch([xsforplot_target fliplr(xsforplot_target)],...
%         [tpos_MED_lower(j,:) fliplr(tpos_MED_upper(j,:))],'g','FaceAlpha',0.5,'EdgeAlpha',0.5); 
    %%%
    
    plot(xsforplot_go,gpos_LOW(j,:),'b'); 
    patch([xsforplot_go fliplr(xsforplot_go)],...
        [gpos_LOW_lower(j,:) fliplr(gpos_LOW_upper(j,:))],'b','FaceAlpha',0.5,'EdgeAlpha',0.5); 

    plot(xsforplot_go,gpos_HIGH(j,:),'r'); 
    patch([xsforplot_go fliplr(xsforplot_go)],...
        [gpos_HIGH_lower(j,:) fliplr(gpos_HIGH_upper(j,:))],'r','FaceAlpha',0.5,'EdgeAlpha',0.5); 

    %%% MEDIUM LIKE
%     plot(xsforplot_go,gpos_MED(j,:),'g'); 
%     patch([xsforplot_go fliplr(xsforplot_go)],...
%         [gpos_MED_lower(j,:) fliplr(gpos_MED_upper(j,:))],'g','FaceAlpha',0.5,'EdgeAlpha',0.5); 
    %%%
    
    title(sprintf('Prior Preference (Prior: %d)',priors{DAYS(j)}.val),'FontSize',16);
    xlabel('Time from Target (ms)','FontSize',14);
    ylabel('Probability','FontSize',14);
    plot([0 0],[-.5 .5],'g--');
    plot([950 950],[-.5 .5],'k--');
    plot([1000 1000],[-.5 .5],'k--');
    ylim([-.5 .5]);
end

%% Prior Preference
figure;
xsforplot_target = -750:100:950;
xsforplot_go = 950 + (50:100:950);
for j = 1:2
    set_of_likes = unique(alldays(prediction_day_indices(j)).tt(:,3));
    
    if length(set_of_likes) > 1
        linds = find(alldays(prediction_day_indices(j)).tt(:,3) == max(set_of_likes));
        hinds = find(alldays(prediction_day_indices(j)).tt(:,3) == min(set_of_likes));
        %minds = find(~ismember(alldays(prediction_day_indices(j)).tt(:,3),[max(set_of_likes) min(set_of_likes)]));
    else
      linds = 1:length(alldays(prediction_day_indices(j)).tt);
      hinds = linds;
      %minds = linds;
    end
        
    for i = 1:size(est_dirs{1}{1},2)
        tpos_LOW(j,i) = nanmean(pri_pref{1}{j}(linds,i));
        [tpos_LOW_lower(j,i),tpos_LOW_upper(j,i)] = boot_bounds(10000,@nanmean,pri_pref{1}{j}(linds,i),2.5,97.5);
        tpos_HIGH(j,i)= nanmean(pri_pref{1}{j}(hinds,i));
        [tpos_HIGH_lower(j,i),tpos_HIGH_upper(j,i)] = boot_bounds(10000,@nanmean,pri_pref{1}{j}(hinds,i),2.5,97.5);
        %%% MEDIUM LIKE STUFF
%         tpos_MED(j,i)= nanmean(pos_like{1}{j}(minds,i));
%         [tpos_MED_lower(j,i),tpos_MED_upper(j,i)] = boot_bounds(10000,@nanmean,pos_like{1}{j}(minds,i),2.5,97.5);
        
    end

    for i = 1:size(est_dirs{2}{1},2)
        gpos_LOW(j,i) = nanmean(pri_pref{2}{j}(linds,i));
        [gpos_LOW_lower(j,i), gpos_LOW_upper(j,i)] = boot_bounds(10000,@nanmean,pri_pref{2}{j}(linds,i),2.5,97.5);
        gpos_HIGH(j,i)= nanmean(pri_pref{2}{j}(hinds,i));
        [gpos_HIGH_lower(j,i), gpos_HIGH_upper(j,i)] = boot_bounds(10000,@nanmean,pri_pref{2}{j}(hinds,i),2.5,97.5);
        %%% MEDIUM LIKE STUFF
%         gpos_MED(j,i)= nanmean(pos_like{2}{j}(minds,i));
%         [gpos_MED_lower(j,i), gpos_MED_upper(j,i)] = boot_bounds(10000,@nanmean,pos_like{2}{j}(minds,i),2.5,97.5);
    end

end

for j = 1:2
    
    subplot(1,2,j); hold on;
    
    plot(xsforplot_target,tpos_LOW(j,:),'b'); 
    patch([xsforplot_target fliplr(xsforplot_target)],...
        [tpos_LOW_lower(j,:) fliplr(tpos_LOW_upper(j,:))],'b','FaceAlpha',0.5,'EdgeAlpha',0.5); 

    plot(xsforplot_target,tpos_HIGH(j,:),'r'); 
    patch([xsforplot_target fliplr(xsforplot_target)],...
        [tpos_HIGH_lower(j,:) fliplr(tpos_HIGH_upper(j,:))],'r','FaceAlpha',0.5,'EdgeAlpha',0.5); 
    
    %%% MEDIUM LIKE
%     plot(xsforplot_target,tpos_MED(j,:),'g'); 
%     patch([xsforplot_target fliplr(xsforplot_target)],...
%         [tpos_MED_lower(j,:) fliplr(tpos_MED_upper(j,:))],'g','FaceAlpha',0.5,'EdgeAlpha',0.5); 
    %%%
    
    plot(xsforplot_go,gpos_LOW(j,:),'b'); 
    patch([xsforplot_go fliplr(xsforplot_go)],...
        [gpos_LOW_lower(j,:) fliplr(gpos_LOW_upper(j,:))],'b','FaceAlpha',0.5,'EdgeAlpha',0.5); 

    plot(xsforplot_go,gpos_HIGH(j,:),'r'); 
    patch([xsforplot_go fliplr(xsforplot_go)],...
        [gpos_HIGH_lower(j,:) fliplr(gpos_HIGH_upper(j,:))],'r','FaceAlpha',0.5,'EdgeAlpha',0.5); 

    %%% MEDIUM LIKE
%     plot(xsforplot_go,gpos_MED(j,:),'g'); 
%     patch([xsforplot_go fliplr(xsforplot_go)],...
%         [gpos_MED_lower(j,:) fliplr(gpos_MED_upper(j,:))],'g','FaceAlpha',0.5,'EdgeAlpha',0.5); 
    %%%
    
    title(sprintf('Prior Preference (Prior: %d)',priors{DAYS(j)}.val),'FontSize',16);
    xlabel('Time from Target (ms)','FontSize',14);
    ylabel('Probability','FontSize',14);
    plot([0 0],[-.5 .5],'g--');
    plot([950 950],[-.5 .5],'k--');
    plot([1000 1000],[-.5 .5],'k--');
    ylim([-.5 .5]);
end

