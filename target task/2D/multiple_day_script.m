monkey = 'MrT';
task = 'UNT2D';
day_list = {'04052013', '04062013'};%,'04072013'}; %{'07102013'};
angle_lims = [-10000 10000];
angle_type = 10;
save_plots = 0;

%% Load files and set up 'session' variable containing all desired days
session = Add_session('MrT',day_list{1},task,0,1,0,[]);
for i = 1:length(day_list)-1
    session = Add_session('MrT',day_list{i+1},'UNT2D',0,1,0,session);
end
orig_place = cd;
%% Use KS_p function to obtain matching neuron lists
[COMPS,~,~,~,ppp] = KS_p(session,0.0025);

alldays = vertcat(session{:}); % Reduce 'session' by concatenating

%% Inputs for PSTH plots
t1_t2 = [-0.80 0]; % Time (s) before and after event
WindowBin = 100; % Size of time bin for PSTH (in ms)
alignment = 'target'; % event to align with ('target','go','reward')

base_day = 1; % Specify which index of day_list will serve as the basis for comparison

%% trial table manipulation
adj_tt = cell(length(day_list),1);
rad_lims = angle_lims*pi/180;
length_tt = zeros(length(day_list),1);
num_feedback_types = zeros(length(day_list),1);
feedback_types = cell(length(day_list),1);
for i = 1:length(day_list)
    
    if rad_lims(1) < rad_lims(2)
        adj_tt{i} = ...
            alldays(i).tt(alldays(i).tt(:,angle_type)>rad_lims(1) & alldays(i).tt(:,angle_type) < rad_lims(2),:);
    else
        adj_tt{i} = ...
            alldays(i).tt(alldays(i).tt(:,angle_type)>rad_lims(1) | alldays(i).tt(:,angle_type) < rad_lims(2),:);
    end
    
    length_tt(i) = length(adj_tt{i});

    feedback_types{i} = unique(adj_tt{i}(:,3));
    num_feedback_types(i) = length(unique(adj_tt{i}(:,3)));
end

large_sets = find(length_tt ~= min(length_tt));
for i = large_sets'   
    if any(length_tt./min(length_tt) > 2)
        rand_trials = randperm(length_tt(i));
        adj_tt{i} = adj_tt{i}(rand_trials(1:num_feedback_types(i)*min(length_tt)),:);
    end
end

% Set color sets
colors{1} = {'k'};
colors{2} = {'r','b'};
colors{3} = {'r','b','g'};

%% Loop & Plot
% Initialize
ax = zeros(length(day_list),1);
fig = cell(length(day_list),1);
s = zeros(length(day_list),1);
rasters = cell(length(day_list),1);
handle = cell(length(day_list),1);
graph_objs = cell(length(day_list),1);
kins = cell(length(day_list),1);

same_neurons = COMPS{base_day}.inds; % Use indices for sorted units
same_chans = COMPS{base_day}.chan;
for i = 1:length(same_chans) 
    same_chans(i,base_day) = alldays(base_day).PMd_units{i}(1);
end
index_pairings = same_neurons;
index_pairings(:,sum(index_pairings,1)==0) = 1:length(index_pairings);

for i = 1:length(same_neurons) % Loop through all neurons of 'base_day'

    neur_pair = same_neurons(i,:); % Find potential unit pairing
    ind_pair = index_pairings(i,:);
    
    if sum(ind_pair==0) == 0 % If a match is found    
    
        for j = 1:length(neur_pair)
            [rasters{j}, handle{j}, graph_objs{j}, kins{j}] = ...
                PSTH_plot(alldays(j).bdfM,alldays(j).PMd_units{ind_pair(j)},...
                          adj_tt{j},t1_t2,WindowBin,alignment,[j-1 0 0]);
            ax(j) = gca;
            fig{j} = get(ax(j),'children');        
        end
        
        % find biggest axis
        graph_days = vertcat(graph_objs{:});
        days_axes = vertcat(graph_days.axis);
        max_axis = find(days_axes(:,end)==max(days_axes(:,end)),1,'first');
        
        new_fig = figure; hold on;
        set(new_fig, 'Position', [100 120 1120 840]);
        
        for j = 1:length(ind_pair)
            s(j) = subplot(2,length(ind_pair),j);
            copyobj(fig{j},s(j));
            
            axis(days_axes(max_axis,:));
            xlabel(graph_objs{j}.xlabel,'FontSize',14);
            ylabel(graph_objs{j}.ylabel,'FontSize',14);
            legend(graph_objs{j}.legend,'Location','Best');
            title(sprintf('(%.1f) Prior (%.1f,%.1f)',same_chans(i,j),...
                graph_objs{j}.prior_mean/pi*180,graph_objs{j}.prior_kappa),'FontSize',16);
            
            close(handle{j});
   
            subplot(2,length(ind_pair),j+length(ind_pair)); hold on;
            
           for k = 1:num_feedback_types(j)
                
                trial_angles_all = alldays(j).tt(alldays(j).tt(:,3)==feedback_types{j}(k),angle_type);
                plot([zeros(length(trial_angles_all),1) cos(trial_angles_all)]',...
                     [zeros(length(trial_angles_all),1) sin(trial_angles_all)]',...
                     'Color',.8*[1 1 1]);
                     %[colors{num_feedback_types(j)}{k} '.-']);
            end
            
            
            trial_angles = cell(num_feedback_types(j),1);
            for k = 1:num_feedback_types(j)
                
                trial_angles{k} = adj_tt{j}(adj_tt{j}(:,3)==feedback_types{j}(k),angle_type);
                plot([zeros(length(trial_angles{k}),1) cos(trial_angles{k})]',...
                     [zeros(length(trial_angles{k}),1) sin(trial_angles{k})]',...
                     'Color',[j-1 0 0]);
                     %[colors{num_feedback_types(j)}{k} '.-']);
            end
            

            
            axis([-1 1 -1 1]);
            %set(gca,'Xtick',[]); set(gca,'Ytick',[]);
            axis off;
            
        end     
        
        if save_plots==1

            cd('C:\Users\limblab\Desktop\figures\4052013_4062013_4072013\go\67p5-112p5\');
            saveas(new_fig,...
                sprintf('PMd_MrT_%s_%d-%d (%0.1f-%0.1f).png',alignment,...
                angle_lims(1),angle_lims(2),same_chans(i,1),same_chans(i,2)),'png');
            close(new_fig);
            cd(orig_place);
            clear new_fig
        end

    end
end

%% Kinematics

time_course = 0:0.750*1000-1;

speed_handle = figure;

for j = 1:length(day_list)
    [rasters{j}, handle{j}, graph_objs{j}, kins{j}] = ...
                PSTH_plot(alldays(j).bdfM,alldays(j).PMd_units{1},...
                          adj_tt{j},[0 0.75],100,'target');
    close(handle{j});
end
                  
for i = 1:length(day_list)  
    
    subplot(1,length(day_list),i); hold on;
    
    if length(graph_objs{i}.legend) > 1
        color_list = {'b','r','g'};
    else
        color_list = {'k'};
    end
    
    for j = 1:length(graph_objs{i}.legend)
    
        av_spd = nanmean(kins{i}{j}.speed,1);
        se_spd = nanstd(kins{i}{j}.speed,1)./sqrt(size(kins{i}{j}.speed,1));
        
        plot(time_course,av_spd , color_list{j}); 
        
        patch([time_course fliplr(time_course)],...
        [av_spd - 1.96*se_spd, fliplr(av_spd + 1.96*se_spd)],...
        color_list{j},'FaceAlpha',0.4,'EdgeAlpha',0.1);
        
    end
end 

%% Get Preferred Directions for tracked neurons
% Specify time bin edges used for tuning curves

%time_bins = [-400 -200 0 200 400 600 800];
tune_range = [700 800];

% Find the tracked units
good_units = find(COMPS{1}.inds(:,2)>0 & COMPS{1}.inds(:,3)>0);
good_units = good_units(2:end);

% Initialize
neurons = cell(length(good_units),length(tune_range)-1);

for i = 1:length(good_units) % Loop through units
     fprintf('%d/%d\n',i,length(good_units));
     for j = 1:length(tune_range)-1 % Loop through time bins
        
        % set time bin edges
        t1 = tune_range(j);
        t2 = tune_range(j+1);
        
        % Run fast_pref_dir.m to calculate tuning curves. 
        [ignore1,cents,tune,tune_se,tune_rast] = ...
            fast_pref_dir(alldays(1).PMd_units,alldays(1).tt(:,10),...
            alldays(1).tt,24,t1,t2,good_units(i));
        
        % Smooth and fix edges on tuning curves
        front_cent = cents(1:3);
        back_cent = cents(end-2:end);
        
        front_tune = repmat(tune(1),1,3); 
        front_se = repmat(tune_se(1),1,3);
        
        back_tune = repmat(tune(end),1,3);
        back_se = repmat(tune_se(end),1,3);
        
        tune_cents = [back_cent-2*pi cents front_cent+2*pi];
        tune_pad = [back_tune';tune;front_tune'];
        se_pad = [back_se';tune_se;front_se'];
        
        interp_cents = tune_cents(1):.01:tune_cents(end);
        tune_interp = interp1(tune_cents,tune_pad,interp_cents);
        se_interp = interp1(tune_cents,se_pad,interp_cents);
        
        smooth_tune = smooth(tune_interp,10);
        smooth_se = smooth(se_interp,10);
        
        wrapped_cents = interp_cents(interp_cents >=0 & interp_cents< 2*pi);
        wrapped_tune = smooth_tune(interp_cents >= 0 & interp_cents <2*pi); 
        wrapped_se = smooth_se(interp_cents >= 0 & interp_cents <2*pi); 
        
        neurons{i,j}.tuning = wrapped_tune;
        neurons{i,j}.tuning_se = wrapped_se;
        neurons{i,j}.raster = tune_rast;
    end
   
end
centers = round(wrapped_cents*1000)./1000;

%% Get Cosine Tunings for tracked neurons
% Specify time bin edges used for tuning curves

%time_bins = [-400 -200 0 200 400 600 800];
tune_range = [700 800];

% Find the tracked units
good_units = find(COMPS{1}.inds(:,2)>0 & COMPS{1}.inds(:,3)>0);
good_units = good_units(2:end);

% Initialize
cos_tunes = cell(length(good_units),1);
max_fr = zeros(length(good_units),1);
min_fr = zeros(length(good_units),1);
normalized = zeros(length(good_units),length(wrapped_cents));
cosfitp = zeros(length(good_units),1);

for i = 1:length(good_units) % Loop through units
    fprintf('%d/%d',i,length(good_units));
        
    % set time bin edges
    t1 = tune_range(1);
    t2 = tune_range(2);

    % Run cosine_tuning.m to calculate tuning curves. 
    [B,fitp,pd,baseline,modul] = cosine_tuning(alldays(1).tt,alldays(1).PMd_units{i},t1,t2,'normal');
    
    theta4tune = (0.01:0.01:2*pi)';
    cosfitp(i) = fitp;
    
    if fitp < 0.05
        % 'log' for poisson, 'identity' for normal
        cos_tunes{i}.fit = glmval(B,[cos(theta4tune) sin(theta4tune)],'identity');
        cos_tunes{i}.norm = (cos_tunes{i}.fit - min(cos_tunes{i}.fit))./(max(cos_tunes{i}.fit)-min(cos_tunes{i}.fit));
        cos_tunes{i}.modul = modul;
        cos_tunes{i}.baseline = baseline;
        fprintf('\n');
    else
        %7 9 12
        cos_tunes{i}.fit = zeros(length(theta4tune),1);
        cos_tunes{i}.norm = zeros(length(theta4tune),1);
        cos_tunes{i}.modul = 0;
        cos_tunes{i}.baseline = 0;
        fprintf(': poorly tuned - dropped\n');
    end
    
    max_fr(i) = max(cos_tunes{i}.fit);
    min_fr(i) = min(cos_tunes{i}.fit);
   
    normalized(i,:) = cos_tunes{i}.norm';
    
end

%% Get Spike Counts
time_bins = -600:100:800;

prediction_day_indices = [2 3];
DAYS_pred = cell(length(prediction_day_indices),1);

for z = 1:length(prediction_day_indices)
    
    prediction_day = prediction_day_indices(z);

    prior_mean = 90; prior_mean = prior_mean/180*pi;
    prior_kappa = 25; 
    day_pred = prediction_day;
    tt_pred = adj_tt{day_pred};

    tracked = COMPS{1}.inds(good_units,day_pred);

    % find rasters for each unit
    unit_pred = cell(length(good_units),1);
    for q = 1:length(tracked)

        [ignore2,Lout,Hout] = ...
            raster_UNT_circ(alldays(day_pred).bdfM,alldays(day_pred).PMd_units{tracked(q)},...
                            tt_pred,time_bins(1)./1000,time_bins(end)./1000);

        Lout(isnan(Lout))=0;
        Hout(isnan(Hout))=0;

        Lrast = zeros(size(Lout,1),length(time_bins)-1);
        Hrast = zeros(size(Hout,1),length(time_bins)-1);
        binsizes = time_bins - time_bins(1); binsizes(1) = 1;
        for v = 1:length(time_bins)-1
            Lrast(:,v) = sum(Lout(:,binsizes(v):binsizes(v+1)),2);
            Hrast(:,v) = sum(Hout(:,binsizes(v):binsizes(v+1)),2);
        end

        linds = find(tt_pred(:,3)==max(tt_pred(:,3)));
        hinds = find(tt_pred(:,3)==min(tt_pred(:,3)));

        full_raster = zeros(length(linds)+length(hinds),length(time_bins)-1);
        full_raster(linds,:) = Lrast;
        full_raster(hinds,:) = Hrast;

        unit_pred{q} = full_raster;
        DAYS_pred{z}{q} = full_raster;

        fprintf('%d/%d\n',q,length(good_units));

    end
end


%% Do population code estimates
DAYS = [2 3];
booted = cell(length(DAYS),1);
for daynumber = 1:length(DAYS)

    prediction_day = DAYS(daynumber);

    %Ps = exp(prior_kappa.*cos(prior_mean-wrapped_cents))./(2.*pi.*besseli(0,prior_kappa));
    Ps = ones(1,length(wrapped_cents));
    %Ps(wrapped_cents > pi) = 0;
    conds = cell(length(time_bins)-1,1);
    conds_cos = cell(length(time_bins)-1,1);
    dt = diff(time_bins(1:2))./1000;

    num_reps = 10000;
    unit_pred = DAYS_pred{prediction_day-1};
    tune_mag = zeros(length(adj_tt{prediction_day}),length(time_bins)-1);
    for bin_num = 1:length(time_bins)-1;

        fprintf('%d/%d\n',bin_num,length(time_bins)-1);

        unit_cell = struct2cell(vertcat(neurons{:,1}));    
        tuning_array = horzcat(unit_cell{1,:})';


        spike_counts = zeros(length(unit_pred),length(adj_tt{prediction_day}));
        for i = 1:length(unit_pred)
            spike_counts(i,:) = unit_pred{i}(:,bin_num)';
        end

        for i = 1:length(adj_tt{prediction_day})
            N = diff(wrapped_cents(1:2));

            calc_post = Ps.*prod(exp(-tuning_array.*dt),1).*prod(tuning_array.^repmat(spike_counts(:,i),1,length(wrapped_cents)),1);

            instant_fr = spike_counts(:,i)/dt;
            multiplier = (instant_fr - min_fr)./(max_fr - min_fr);
            multiplier(isnan(multiplier)) = 0;
            multiplier(isinf(multiplier)) = 0;


            mult_mat = repmat(multiplier,1,size(normalized,2));

            conds{bin_num}(i,:) = calc_post./(sum(calc_post)*N);
            conds_cos{bin_num}(i,:) = sum(mult_mat.*normalized);

            tune_mag(i,bin_num) = max(conds_cos{bin_num}(i,:),[],2) - min(conds_cos{bin_num}(i,:),[],2);

        end     
    end

    booted{daynumber} = bootstrp(10000,@mean,tune_mag);

end
%%

conds = C3;
dayind = 3;
figure; hold on;

end_peak = zeros(length(conds),size(conds{1},1));
targ_peak = zeros(length(conds),size(conds{1},1));
v_estimate = zeros(length(conds),size(conds{1},1));
pop_est = zeros(length(conds),size(conds{1},1));
end_like = zeros(length(conds),size(conds{1},1));
targ_like = zeros(length(conds),size(conds{1},1));
popvar = zeros(length(conds),size(conds{1},1));
popmean = zeros(length(conds),size(conds{1},1));
popconc = zeros(length(conds),size(conds{1},1));

for i = 1:length(conds)
    
    fprintf('%d/%d\n',i,length(conds));
    
    targ_shift = zeros(size(conds{i}));
    end_shift = zeros(size(conds{i}));
    
    for j = 1:size(conds{i},1)
        
        targ_pos = adj_tt{dayind}(j,9);
        end_pos = adj_tt{dayind}(j,10);
        
        targ_ind = find(wrapped_cents > targ_pos,1,'first');
        end_ind = find(wrapped_cents > end_pos,1,'first');
        
        targ_shift(j,:) = [conds{i}(j,targ_ind:end) conds{i}(j,1:targ_ind-1)];
        end_shift(j,:) = [conds{i}(j,end_ind:end) conds{i}(j,1:end_ind-1)];
        
        if ~isnan(max(conds{i}(j,:)))
            pop_est(i,j) = wrapped_cents(conds{i}(j,:)==max(conds{i}(j,:)));
            end_like(i,j) = end_shift(j,1);
            targ_like(i,j) = targ_shift(j,1);
            
            [samps] = sample_from_pdf(wrapped_cents,conds{i}(j,:),1000);
            %sum(isnan(samps))
            popvar(i,j) = circ_var(samps);
            popmean(i,j) = circ_mean(samps);
            [g1,g2,popconc(i,j)] = calc_disp(samps);
            
        else
            pop_est(i,j) = nan;
            end_like(i,j) = nan;
            targ_like(i,j) = nan;
        end
    
    end
    
    thetas = wrapped_cents; 
    thetas(thetas > pi) = thetas(thetas > pi)-2*pi;

    targ_aligned = [targ_shift(:, thetas < 0) targ_shift(:,thetas > 0)]; 
    end_aligned = [end_shift(:, thetas < 0) end_shift(:,thetas > 0)]; 
    thetas_ordered = [thetas(thetas < 0) thetas(thetas > 0)];

%     figure; plot(thetas_ordered, nanmean(targ_aligned,1),'r'); hold on ;
%     plot(thetas_ordered, nanmean(end_aligned,1),'b');
    subplot(2,length(conds)./2,i);
    %plot(thetas_ordered, nanmean(targ_aligned),'b');
    %plot(thetas_ordered, nanmean(end_aligned),'b');
    plot(thetas_ordered, nanmean(end_aligned(adj_tt{dayind}(:,3)==min(adj_tt{dayind}(:,3)),:),1),'r');
    hold on;
    plot(thetas_ordered, nanmean(end_aligned(adj_tt{dayind}(:,3)==max(adj_tt{dayind}(:,3)),:),1),'b');

    title(sprintf('bin: %d',i));

    for q = 1:length(end_aligned)
        if sum(isnan(end_aligned(q,:))) == 0
            end_peak(i,q) = thetas_ordered(end_aligned(q,:)==max(end_aligned(q,:)));
            targ_peak(i,q) = thetas_ordered(targ_aligned(q,:)==max(targ_aligned(q,:)));
            mu_v = sum(thetas_ordered.*N.*end_aligned(q,:));
            v_estimate(i,q) = 1/max(end_aligned(q,:));%sum((thetas_ordered - mu_v).^2.*end_aligned(q,:));
            
        else
            end_peak(i,q) = nan;
            targ_peak(i,q) = nan;
            v_estimate(i,q) = nan;
        end
    end
    %clc;
end

%%
tune_mag = zeros(size(C2C{1},1),length(C2C));
for i = 1:length(C2C)
    tune_mag(:,i) = max(C2C{i},[],2) - min(C2C{i},[],2);
end

tune_mag2 = zeros(size(C3C{1},1),length(C3C));
for i = 1:length(C3C)
    tune_mag2(:,i) = max(C3C{i},[],2) - min(C3C{i},[],2);
end

figure; plot(nanmean(tune_mag,1),'b');
hold on; plot(nanmean(tune_mag2,1),'g');
% figure; hold on;
% plot(nanmean(tune_mag(linds,:),1),'b');
% plot(nanmean(tune_mag(hinds,:),1),'r');


%% Do tuning probabilities

prediction_day = 2;

prior_mean = 90; prior_mean = prior_mean/180*pi;
prior_kappa = 7; 

%Ps = exp(prior_kappa.*cos(prior_mean-wrapped_cents))./(2.*pi.*besseli(0,prior_kappa));
Ps = ones(1,length(wrapped_cents));
%Ps(wrapped_cents > pi) = 0;
conds = cell(length(time_bins)-1,1);
dt = diff(time_bins(1:2))./1000;

unit_pred = DAYS_pred{prediction_day-1};

for bin_num = 1:length(time_bins)-1;
    unit_cell = struct2cell(vertcat(neurons{:,bin_num}));    
    tuning_array = horzcat(unit_cell{1,:})';
    
    spike_counts = zeros(length(unit_pred),length(adj_tt{prediction_day}));
    for i = 1:length(unit_pred)
        spike_counts(i,:) = unit_pred{i}(:,bin_num)';
    end
    
    for i = 1:length(adj_tt{prediction_day})
        N = diff(wrapped_cents(1:2));
        
        calc_post = Ps.*prod(exp(-tuning_array.*dt),1).*prod(tuning_array.^repmat(spike_counts(:,i),1,length(wrapped_cents)),1);

        conds{bin_num}(i,:) = calc_post./(sum(calc_post)*N);
        
    end     
end


    