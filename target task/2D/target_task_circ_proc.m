%-% Setup/User Input
computer = 'work';
monkey = 'Mihili';
MO = '07';  %%7-19-2013
DA = '19';
YE = '2013';
TASK = 'UNT2D';
fileindx = '001-01';

sortedM1 = 1;
Behavior_plots = 0;
behavior_only = 0; % 1 if only loading behavior

%- BEHAVIOR PARAMETERS ----------------------------------------------------
shift_m = 5.5; % Shift mean
shift_v = 20; % Shift variance
numslices = 5;

yoff = -3; % Global y offset of center target
target_rad = 10; % Target radius
target_size = 2; % Target size

% Screen offsets
x_off = -3; 
y_off = 33;

%-% Load files
fprintf('Loading files...\n');
modaye = [MO DA YE];
orig_place = cd;

% Go to folder with BDFs %
if strcmp(computer,'work')
    cd(sprintf('C:\\Users\\limblab\\Desktop\\%s_bdf',monkey));
else   
    cd('../../');
    cd([cd '/bdf'])
end

% Load M1/Behavior bdf
if sortedM1==1
    M1_file = sprintf('bdf_sorted_%s_M1_%s_%s_%s.mat',monkey,modaye,TASK,fileindx);
    load(M1_file);
    bdfM = bdf; clear bdf;
else
    M1_file = sprintf('bdf_%s_M1_%s_%s_%s.mat',monkey,modaye,TASK,fileindx(1:3));
    load(M1_file);
    bdfM = bdf; clear bdf;
end

% Load PMd bdf
if behavior_only ~= 1
    PMd_file = sprintf('bdf_sorted_%s_PMd_%s_%s_%s.mat',monkey,modaye,TASK,fileindx);
    load(PMd_file);
    bdfP = bdf; clear bdf;
end

% Go back to original folder and clean up
cd(orig_place);
clear M1_file PMd_file fileindx M1_units PMd_units
%-% Get Unit information
id_func = @(x) x(1);
if behavior_only ~= 1
    if sortedM1==1
        M1_units = spiketrains(bdfM,1);
        M1_ids = cellfun(id_func,M1_units);
    else
        M1_units = [];
        M1_ids = [];
    end
    PMd_units = spiketrains(bdfP,1);
    PMd_ids = cellfun(id_func,PMd_units);
end
%-% Get Behavior
fprintf('Getting behavior information...\n');
%[tt] = getTT_UNT_circ_OLD(bdfM);
[tt,tts,priors] = getTT_UNT_circ(bdfM);
[BDF ,comp_tt, tt_labels, Linds, Hinds, kin, slices] = UNT_circ_behavior(bdfM,tt,x_off,y_off,numslices);

c_tt = cell(length(tts),1);
slice_table = cell(length(tts),1);
for i = 1:length(tts)
    [~,c_tt{i},~,~,~,~,slice_table{i}] = UNT_circ_behavior(bdfM,tts{i},x_off,y_off,numslices);
end

fprintf('Correcting for photo-delay and catch trials...\n');
%-% Photo Correction
orig_comp_tt = comp_tt; 
[corrected_ts,lags] = photo_correct(comp_tt,bdfM,[5 6 7],'mean');
comp_tt(:,[5 6 7]) = corrected_ts;

if exist('c_tt','var')
    
    o_c_tt = cell(length(c_tt),1);
    c_ts = cell(length(c_tt),1);
    for i = 1:length(c_tt)
        
        o_c_tt{i} = c_tt{i};
        [c_ts{i}] = photo_correct(c_tt{i},bdfM,[5 6 7],'mean');
        c_tt{i}(:,[5 6 7]) = c_ts{i};
        
    end
end

%-% Drop Catch trials
[nocatch_tt, catch_tt] = dropCatchTrials(comp_tt,0.5);
comp_tt = nocatch_tt; 
good_inds = cell(length(tts),1);

if exist('c_tt','var')
    
    good_tt = cell(length(c_tt),1);
    bad_tt = cell(length(c_tt),1);
    for i = 1:length(c_tt)
       
       [good_tt{i},bad_tt{i},good_inds{i}] = dropCatchTrials(c_tt{i},0.5);
       c_tt{i} = good_tt{i};
       slice_table{i} = slice_table{i}(good_inds{i},:);
       
   end
end

%-% Set up new variables for cross-day comparison
if ~exist('session_nokin','var'); session_nokin = cell(12,31); end
if ~exist('bdfP','var'); bdfP = 0; end
if ~exist('PMd_units','var'); PMd_units = 0; end

session_nokin{str2double(MO),str2double(DA)}.tt = comp_tt;
session_nokin{str2double(MO),str2double(DA)}.PMd_units = PMd_units;
session_nokin{str2double(MO),str2double(DA)}.M1_units = M1_units;

if exist('c_tt','var')
    
    clear session_nokin;
    session_nokin = cell(length(c_tt),1);
    for i = 1:length(c_tt)
        session_nokin{i}.tt = c_tt{i};
        session_nokin{i}.PMd_units = PMd_units;
        session_nokin{i}.M1_units = M1_units;
        session_nokin{i}.slices = slice_table{i};
    end
end
alldays_nokin = vertcat(session_nokin{:}); % Reduce 'session' by concatenating

%-% Set up no-kin variables for cross-day comparison
if ~exist('session','var'); session = cell(12,31); end
if ~exist('bdfP','var'); bdfP = 0; end
if ~exist('PMd_units','var'); PMd_units = 0; end

session{str2double(MO),str2double(DA)}.bdfM = bdfM;
session{str2double(MO),str2double(DA)}.bdfP = bdfP;
session{str2double(MO),str2double(DA)}.tt = comp_tt;
session{str2double(MO),str2double(DA)}.PMd_units = PMd_units;
session{str2double(MO),str2double(DA)}.M1_units = M1_units;
session{str2double(MO),str2double(DA)}.kin = BDF;
session{str2double(MO),str2double(DA)}.vars.Linds = Linds;
session{str2double(MO),str2double(DA)}.vars.Hinds = Hinds;
session{str2double(MO),str2double(DA)}.slices = slices;

if exist('c_tt','var')
    
    clear session;
    session = cell(length(c_tt),1);
    for i = 1:length(c_tt)
        session{i}.bdfM = bdfM;
        session{i}.bdfP = bdfP;
        session{i}.tt = c_tt{i};
        session{i}.PMd_units = PMd_units;
        session{i}.M1_units = M1_units;
        session{i}.kin = BDF;
        session{i}.vars.Linds = Linds;
        session{i}.vars.Hinds = Hinds;
        session{i}.slices = slice_table{i};
    end
end

alldays = vertcat(session{:}); % Reduce 'session' by concatenating

fprintf('Done\n');

%% Do Plots 
if Behavior_plots == 1 && ~exist('c_tt','var')
    plot_tt = comp_tt;
    Behavior_Plots_unc_circ;
elseif Behavior_plots == 1 && exist('c_tt','var')
    
    for i = 1:length(c_tt)
        plot_tt = c_tt{i};
        Behavior_Plots_unc_circ;
    end
end
%% Add waveform information
bdf_in_use = bdfM; brain_reg = 'M1';
bdf_add_waveform; bdfM = bdf_in_use;
if isfield(bdfM.units,'wave'); M1_waves = wave_shapes(bdfM); end;
bdf_in_use = bdfP; brain_reg = 'PMd';
bdf_add_waveform; bdfP = bdf_in_use;
if isfield(bdfP.units,'wave'); PMd_waves = wave_shapes(bdfP); end;
%% Loop & Plot Raster M1
t1_t2 = [-0.8 .8]; % Time (s) before and after event
alignment = 'target'; % event to align with ('target','go','reward')

save_plotsM1 = 1;
day_spikes = alldays(1).M1_units;
day_tt = alldays(1).tt;

mkdir(sprintf('C:\\Users\\limblab\\Desktop\\figures\\%s\\M1\\%s\\%s\\prior\\',monkey,alignment,modaye));
if sortedM1

    %mkdir(sprintf('C:\\Users\\limblab\\Desktop\\figures\\%s\\M1\\%s\\%s\\',monkey,alignment,modaye));
    for i = 1:length(day_spikes)
      
        chanind = ['C' num2str(floor(day_spikes{i}(1))) 'U',...
            num2str(10*(day_spikes{i}(1)-floor(day_spikes{i}(1))))];
        
        fprintf('%s\n',chanind);
        
        [LRast,HRast,fig_hand] = raster_plot(day_spikes{i},day_tt,t1_t2,...
            alignment,1,'endpoint',sprintf('Unit %s',chanind));
   
        if save_plotsM1==1
            cd(sprintf('C:\\Users\\limblab\\Desktop\\figures\\%s\\M1\\%s\\%s\\prior\\',monkey,alignment,modaye));
                saveas(fig_hand,...
                    sprintf('%s_%s_%s_%s_%s.png',monkey,TASK,modaye,alignment,chanind),'png');
                close(fig_hand);
                cd(orig_place);
                clear fig_hand
        end
    
    end
end

%% Loop & Plot Raster PMd
t1_t2 = [0 .6]; % Time (s) before and after event
alignment = 'go'; % event to align with ('target','go','reward')

save_plotsPMd = 0;
day_spikes = alldays(1).PMd_units;
day_tt = alldays(2).tt; %alldays(1).tt;

mkdir(sprintf('C:\\Users\\limblab\\Desktop\\figures\\%s\\PMd\\%s\\%s\\',monkey,alignment,modaye));
psthPMD = zeros(length(day_spikes),9);
for i = 104%1:length(day_spikes)
    
    chanind = ['C' num2str(floor(day_spikes{i}(1))) 'U',...
        num2str(10*(day_spikes{i}(1)-floor(day_spikes{i}(1))))];

    fprintf('%s\n',chanind);

    [Rast,allinds,fig_hand] = raster_plot(day_spikes{i},day_tt,t1_t2,...
        alignment,1,'endpoint',sprintf('Unit %s',chanind));
    
%     [L_Rast,allinds,fig_hand] = raster_plot(day_spikes{i},day_tt,t1_t2,...
%         alignment,1,'trials',sprintf('Unit %s',chanind));

    if save_plotsPMd==1
        cd(sprintf('C:\\Users\\limblab\\Desktop\\figures\\%s\\PMd\\%s\\%s\\',monkey,alignment,modaye));
            saveas(fig_hand,...
                sprintf('%s_%s_%s_%s_%s.png',monkey,TASK,modaye,alignment,chanind),'png');
            close(fig_hand);
            cd(orig_place);
            clear fig_hand
    end
    
%     PSTH_HI(:,:,i) = bin_array(Rast{1},size(Rast{1},1),30);
%     PSTH_LO(:,:,i) = bin_array(Rast{2},size(Rast{2},1),30);
%     
%     PSTH_HI(isnan(PSTH_HI)) = 0;
%     PSTH_LO(isnan(PSTH_LO)) = 0;
%     LowRast = bin_array(Rast{1},1,9); HighRast = bin_array(Rast{2},1,9);
%     figure; plot(LowRast,'b'); hold on; plot(HighRast,'r'); pause; cla;
%     psthPMD(i,:) = bin_array(vertcat(Rast{:}),1,9);
    
end

%% Get Cosine Tunings for tracked neurons
% Specify time bin edges used for tuning curves

tune_range = [0 250];
tune_align = 'go';

co_units = alldays(1).PMd_units;

% Initialize
cos_tunes = cell(length(co_units),1);
max_fr = zeros(length(co_units),1);
min_fr = zeros(length(co_units),1);
normalized = zeros(length(co_units),length(wrapped_cents));
cosfitp = zeros(length(co_units),1);

good_neuron_count = 0;
for i = 1:length(co_units) % Loop through units
    fprintf('%d/%d',i,length(co_units));
        
    % set time bin edges
    t1 = tune_range(1);
    t2 = tune_range(2);

    % Run cosine_tuning.m to calculate tuning curves. 
    [B,fitp,pd,baseline,modul] = cosine_tuning(alldays(1).tt,co_units{i},t1,t2,'normal',tune_align);
    
    theta4tune = (0.01:0.01:2*pi)';
    cosfitp(i) = fitp;
    
    if fitp < 0.05
        % 'log' for poisson, 'identity' for normal
        cos_tunes{i}.fit = glmval(B,[cos(theta4tune) sin(theta4tune)],'identity');
        cos_tunes{i}.norm = (cos_tunes{i}.fit - min(cos_tunes{i}.fit))./(max(cos_tunes{i}.fit)-min(cos_tunes{i}.fit));
        cos_tunes{i}.modul = modul;
        cos_tunes{i}.baseline = baseline;
        fprintf('\n');
        good_neuron_count = good_neuron_count + 1;
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
%     pause; close;
end
fprintf('\n%d/%d units used\n',good_neuron_count,length(co_units));

%% Get Preferred Directions and non-parameterized tuning curves for tracked neurons
% Specify time bin edges used for tuning curves

tune_range = [0 250];
tune_align = 'go';
brain_area = 'PMd';

% Initialize
co_units = eval(sprintf('alldays(1).%s_units',brain_area));
if strcmp(brain_area,'M1'); co_units(end-1:end) =[]; end
neurons = cell(length(co_units),length(tune_range)-1);

for i = 1:length(co_units) % Loop through units
     fprintf('%d/%d\n',i,length(co_units));
     
    % set time bin edges
    t1 = tune_range(1);
    t2 = tune_range(2);

    % Run fast_pref_dir.m to calculate tuning curves. 
%     [ignore1,cents,tune,tune_se,tune_rast] = ...
%         fast_pref_dir(co_units,alldays(1).tt(:,10),...
%         alldays(1).tt,16,t1,t2,tune_align,i);

    [cents,tune,tune_se] = co_tuning(co_units,alldays(1).tt(:,10),alldays(1).tt,t1,t2,tune_align,i);
    
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

    smooth_tune = smooth(tune_interp,100);
    smooth_se = smooth(se_interp,100);

    wrapped_cents = interp_cents(interp_cents >=0 & interp_cents< 2*pi);
    wrapped_tune = smooth_tune(interp_cents >= 0 & interp_cents <2*pi); 
    wrapped_se = smooth_se(interp_cents >= 0 & interp_cents <2*pi); 

    neurons{i}.tuning = wrapped_tune;
    neurons{i}.tuning_se = wrapped_se;
    %neurons{i}.raster = tune_rast;
   
end
centers = round(wrapped_cents*1000)./1000;

%% Get Spike Counts
time_bins = -500:150:1000;
time_align = 'target';
brain_area = 'PMd';
prior_mean = 90;

prediction_day_indices = [2 3];
DAYS_pred = cell(length(prediction_day_indices),1);
prior_mean = prior_mean/180*pi;
for z = 1:length(prediction_day_indices)
   
    prediction_day = prediction_day_indices(z);
    day_units = eval(sprintf('alldays(%d).%s_units',prediction_day,brain_area));
    if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end

    day_pred = prediction_day;
    tt_pred = alldays(day_pred).tt;

    % find rasters for each unit
    unit_pred = cell(length(day_units),1);
    for q = 1:length(day_units)
     
        [Lout,Hout] = raster_plot(day_units{q},tt_pred,[time_bins(1)./1000 time_bins(end)./1000],...
            time_align,0,'none');
        
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

        fprintf('%d/%d\n',q,length(day_units));

    end
end

%% Do population code estimates
DAYS = [2, 3];
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
for daynumber = 1:length(DAYS)

    prediction_day = DAYS(daynumber);

    %Ps = exp(prior_kappa.*cos(prior_mean-wrapped_cents))./(2.*pi.*besseli(0,prior_kappa));
    Ps = ones(1,length(wrapped_cents));
    %Ps(wrapped_cents > pi) = 0;
    conds{daynumber} = cell(length(time_bins)-1,1);
    cond_sum{daynumber} = cell(length(time_bins)-1,1);
    conds_cos{daynumber} = cell(length(time_bins)-1,1);
    estimated_dir{daynumber} = zeros(length(alldays(prediction_day).tt),length(time_bins)-1);
    estimate_error_dir{daynumber} = zeros(length(alldays(prediction_day).tt),length(time_bins)-1);
    estimate_error_cen{daynumber} = zeros(length(alldays(prediction_day).tt),length(time_bins)-1);
    estimate_conf{daynumber} = zeros(length(alldays(prediction_day).tt),length(time_bins)-1);
    actual_dir{daynumber} = zeros(length(alldays(prediction_day).tt),1);
    actual_cen{daynumber} = zeros(length(alldays(prediction_day).tt),1);
    
    dt = diff(time_bins(1:2))./1000;

    num_reps = 10000;
    unit_pred = DAYS_pred{daynumber};
    tune_mag = zeros(length(alldays(prediction_day).tt),length(time_bins)-1);
    for bin_num = 1:length(time_bins)-1;

        fprintf('day %d: %d/%d\n',daynumber,bin_num,length(time_bins)-1);

        unit_cell = struct2cell(vertcat(neurons{:,1}));    
        tuning_array = horzcat(unit_cell{1,:})';

        spike_counts = zeros(length(unit_pred),length(alldays(prediction_day).tt));
        for i = 1:length(unit_pred)
            spike_counts(i,:) = unit_pred{i}(:,bin_num)';
        end

        for i = 1:length(alldays(prediction_day).tt)
            N = diff(wrapped_cents(1:2));

%             tuning_array = tuning_array(1:20,:);
%             spike_counts = spike_counts(1:20,:);
            
            paren = Ps.*prod(exp(-tuning_array.*dt),1);
            ta_pow = tuning_array.^repmat(spike_counts(:,i),1,length(wrapped_cents));
            ta_pow = ta_pow./repmat(max(ta_pow,[],2),1,size(ta_pow,2));
            ta_pow(isnan(ta_pow))=1;
            calc_post = paren.*prod(ta_pow,1);
            %calc_post = Ps.*prod(exp(-tuning_array.*dt),1).*prod(tuning_array.^repmat(spike_counts(:,i),1,length(wrapped_cents)),1);
%             paren = repmat(spike_counts(:,i),1,length(wrapped_cents)).*log(tuning_array) - tuning_array.*dt;
%             paren(isinf(paren))=0; paren(isnan(paren))=0;
%             exponent = nansum(paren,1);
%             calc_post = exp(exponent);
%             
%             USE FOR COSINE TUNED PREDICTIONS
%             instant_fr = spike_counts(:,i)/dt;
%             multiplier = (instant_fr - min_fr)./(max_fr - min_fr);
%             multiplier(isnan(multiplier)) = 0;
%             multiplier(isinf(multiplier)) = 0;
% 
%             mult_mat = repmat(multiplier,1,size(normalized,2));
%             conds_cos{daynumber}{bin_num}(i,:) = sum(mult_mat.*normalized);
%             tune_mag(i,bin_num) = max(conds_cos{daynumber}{bin_num}(i,:),[],2) - min(conds_cos{daynumber}{bin_num}(i,:),[],2);   
            
            conds{daynumber}{bin_num}(i,:) = calc_post./(sum(calc_post)*N);
           
            with_nans = conds{daynumber}{bin_num}(i,:);
            with_nans(isinf(with_nans))=0;
            with_nans(isnan(with_nans))=0;
            samps = sample_from_pdf(wrapped_cents,with_nans,10000);
            
            %estimated_dir{daynumber}(i,bin_num) = mean_angle(samps,'rads');
            estimated_dir{daynumber}(i,bin_num) = wrapped_cents(find(with_nans == max(with_nans),1,'first'));
            
            actual_dir{daynumber}(i) = alldays(prediction_day).tt(i,10);
            actual_cen{daynumber}(i) = alldays(prediction_day).tt(i,9);
            
            estimate_error_dir{daynumber}(i,bin_num) = ...
                angle_diff(estimated_dir{daynumber}(i,bin_num),actual_dir{daynumber}(i));
            
            estimate_error_cen{daynumber}(i,bin_num) = ...
                angle_diff(estimated_dir{daynumber}(i,bin_num),actual_cen{daynumber}(i));
            
            estimate_conf{daynumber}(i,bin_num) = max(with_nans);
        end     
    end

    booted{daynumber} = bootstrp(10000,@mean,tune_mag);

end
%%
%Booted errorbar.
times = (time_bins(2:end)+time_bins(1:end-1))./2;
m = cell(length(booted),1);
color_pal = {'r','b','g'};
figure; hold on;

pri_struct = vertcat(priors{:});
pri_list = horzcat(pri_struct(DAYS).val);
leg_struct = cell(length(booted),1);

for i = 1:length(booted)
    m{i} = mean(booted{i},1);
    plot(times,m{i},color_pal{i});
end

for i = 1:length(booted)

    e_mat = sort(booted{i});

    l_e_ind = round(size(booted{i},1)*0.025);
    h_e_ind = size(booted{i},1) - l_e_ind + 1;


    patch([times fliplr(times)],[e_mat(l_e_ind,:) fliplr(e_mat(h_e_ind,:))],color_pal{i},'FaceAlpha',0.5,'EdgeAlpha',0.5);
    
    leg_struct{i} = ['Prior: ' num2str(pri_list(i))];
end

legend(leg_struct);

%% Plot Error
figure; hold on;
xsforplot = 0.5*(time_bins(1:end-1)+time_bins(2:end));
for i = 1:length(DAYS);
    
    linds = find(alldays(DAYS(i)).tt(:,3)==max(alldays(DAYS(i)).tt(:,3)));
    hinds = find(alldays(DAYS(i)).tt(:,3)==min(alldays(DAYS(i)).tt(:,3)));
    
    [dml, dul_l, dll_l] = circ_mean(abs(estimate_error_dir{i}(linds,:)));
    [dmh, dul_h, dll_h] = circ_mean(abs(estimate_error_dir{i}(hinds,:)));
    
    [cml, cul_l, cll_l] = circ_mean(abs(estimate_error_cen{i}(linds,:)));
    [cmh, cul_h, cll_h] = circ_mean(abs(estimate_error_cen{i}(hinds,:)));
    
    subplot(1,length(DAYS),i); hold on;
    patch([xsforplot fliplr(xsforplot)],[dll_l fliplr(dul_l)],'b','FaceAlpha',0.5,'EdgeAlpha',0.5);
    patch([xsforplot fliplr(xsforplot)],[dll_h fliplr(dul_h)],'r','FaceAlpha',0.5,'EdgeAlpha',0.5);
    plot(xsforplot,dml,'b');
    plot(xsforplot,dmh,'r');
    title(sprintf('End Position Error (Prior: %d)',priors{DAYS(i)}.val),'FontSize',16);
    xlabel('Time from Target (ms)','FontSize',14);
    ylabel('abs(prediction error)','FontSize',14);
    plot([0 0],[0 pi],'g--');
    ylim([0 pi]);
    
%     subplot(2,length(DAYS),i+length(DAYS)); hold on;
%     patch([xsforplot fliplr(xsforplot)],[cll_l fliplr(cul_l)],'b','FaceAlpha',0.5,'EdgeAlpha',0.5);
%     patch([xsforplot fliplr(xsforplot)],[cll_h fliplr(cul_h)],'r','FaceAlpha',0.5,'EdgeAlpha',0.5);
%     plot(xsforplot,cml,'b');
%     plot(xsforplot,cmh,'r');
%     title(sprintf('Centroid Position Error (Prior: %d)',priors{DAYS(i)}.val),'FontSize',16);
%     ylim([0 pi]);
    
end

%% Plot Confidence for two epochs
figure; hold on;
xsforplot_target = -750:100:950;
xsforplot_go = 950 + (50:100:950);
for i = 1:length(DAYS);
    
    linds = find(alldays(DAYS(i)).tt(:,3)==max(alldays(DAYS(i)).tt(:,3)));
    hinds = find(alldays(DAYS(i)).tt(:,3)==min(alldays(DAYS(i)).tt(:,3)));
    
    % TARGET epoch  
    conl.t = mean(conf_targ{i}(linds,:));
    conl_u.t = conl.t + 1.95*std(conf_targ{i}(linds,:))./sqrt(size(conf_targ{i}(linds,:),1));
    conl_l.t = conl.t - 1.95*std(conf_targ{i}(linds,:))./sqrt(size(conf_targ{i}(linds,:),1));
    
    conh.t = mean(conf_targ{i}(hinds,:));
    conh_u.t = conh.t + 1.95*std(conf_targ{i}(hinds,:))./sqrt(size(conf_targ{i}(hinds,:),1));
    conh_l.t = conh.t - 1.95*std(conf_targ{i}(hinds,:))./sqrt(size(conf_targ{i}(hinds,:),1));
    
    % GO epoch
    conl.g = mean(conf_go{i}(linds,:));
    conl_u.g = conl.g + 1.95*std(conf_go{i}(linds,:))./sqrt(size(conf_go{i}(linds,:),1));
    conl_l.g = conl.g - 1.95*std(conf_go{i}(linds,:))./sqrt(size(conf_go{i}(linds,:),1));

    conh.g = mean(conf_go{i}(hinds,:));
    conh_u.g = conh.g + 1.95*std(conf_go{i}(hinds,:))./sqrt(size(conf_go{i}(hinds,:),1));
    conh_l.g = conh.g - 1.95*std(conf_go{i}(hinds,:))./sqrt(size(conf_go{i}(hinds,:),1));
    
    %PLOT Target
    subplot(1,length(DAYS),i); hold on;
    patch([xsforplot_target fliplr(xsforplot_target)],[conl_l.t fliplr(conl_u.t)],'b','FaceAlpha',0.5,'EdgeAlpha',0.5);
    patch([xsforplot_target fliplr(xsforplot_target)],[conh_l.t fliplr(conh_u.t)],'r','FaceAlpha',0.5,'EdgeAlpha',0.5);
    plot(xsforplot_target,conl.t,'b');
    plot(xsforplot_target,conh.t,'r');
    
    %PLOT Go
    patch([xsforplot_go fliplr(xsforplot_go)],[conl_l.g fliplr(conl_u.g)],'b','FaceAlpha',0.5,'EdgeAlpha',0.5);
    patch([xsforplot_go fliplr(xsforplot_go)],[conh_l.g fliplr(conh_u.g)],'r','FaceAlpha',0.5,'EdgeAlpha',0.5);
    plot(xsforplot_go,conl.g,'b');
    plot(xsforplot_go,conh.g,'r');
    plot([950 950],[0 10],'k--');
    plot([1000 1000],[0 10],'k--');
    
    title(sprintf('Prediction Confidence (Prior: %d)',priors{DAYS(i)}.val),'FontSize',16);
    xlabel('Time from Target (ms)','FontSize',14);
    ylabel('confidence','FontSize',14);
    plot([0 0],[0 10],'g--');
    ylim([0 10]);

end


%% Plot Absolute Error for two epochs
figure; hold on;
xsforplot_target = -750:100:950;
xsforplot_go = 950 + (50:100:950);
for i = 1:length(DAYS);
    
    linds = find(alldays(DAYS(i)).tt(:,3)==max(alldays(DAYS(i)).tt(:,3)));
    hinds = find(alldays(DAYS(i)).tt(:,3)==min(alldays(DAYS(i)).tt(:,3)));
    
    % TARGET epoch      
    
    [dml.t, dul_l.t, dll_l.t] = circ_mean(abs(eed_target{i}(linds,:)));
    targ_L_std = circ_std(abs(eed_target{i}(linds,:)));
    [dmh.t, dul_h.t, dll_h.t] = circ_mean(abs(eed_target{i}(hinds,:)));
    targ_H_std = circ_std(abs(eed_target{i}(hinds,:)));
    
    % GO epoch
    [dml.g, dul_l.g, dll_l.g] = circ_mean(abs(eed_go{i}(linds,:)));
    go_L_std = circ_std(abs(eed_go{i}(linds,:)));
    [dmh.g, dul_h.g, dll_h.g] = circ_mean(abs(eed_go{i}(hinds,:)));
    go_H_std = circ_std(abs(eed_go{i}(hinds,:)));
    
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
    
    linds = find(alldays(DAYS(i)).tt(:,3)==max(alldays(DAYS(i)).tt(:,3)));
    hinds = find(alldays(DAYS(i)).tt(:,3)==min(alldays(DAYS(i)).tt(:,3)));
    
    % TARGET epoch  
    [cml.t, cul_l.t, cll_l.t] = circ_mean(eec_target{i}(linds,:));
    [cmh.t, cul_h.t, cll_h.t] = circ_mean(eec_target{i}(hinds,:));
   
    [dml.t, dul_l.t, dll_l.t] = circ_mean(eed_target{i}(linds,:));
    [dmh.t, dul_h.t, dll_h.t] = circ_mean(eed_target{i}(hinds,:));
    
%     [dml.t, dul_l.t, dll_l.t] = circ_mean(eed_target{i}(linds,:));
%     [dmh.t, dul_h.t, dll_h.t] = circ_mean(eed_target{i}(hinds,:));
%     
%     [cml.t, cul_l.t, cll_l.t] = circ_mean(eec_target{i}(linds,:));
%     [cmh.t, cul_h.t, cll_h.t] = circ_mean(eec_target{i}(hinds,:));    
    
    % GO epoch
    [dml.g, dul_l.g, dll_l.g] = circ_mean(eed_go{i}(linds,:));
    [dmh.g, dul_h.g, dll_h.g] = circ_mean(eed_go{i}(hinds,:));
    
    [cml.g, cul_l.g, cll_l.g] = circ_mean(eec_go{i}(linds,:));
    [cmh.g, cul_h.g, cll_h.g] = circ_mean(eec_go{i}(hinds,:));
    
%     [dml.g, dul_l.g, dll_l.g] = circ_mean(eed_go{i}(linds,:));
%     [dmh.g, dul_h.g, dll_h.g] = circ_mean(eed_go{i}(hinds,:));
%     
%     [cml.g, cul_l.g, cll_l.g] = circ_mean(eec_go{i}(linds,:));
%     [cmh.g, cul_h.g, cll_h.g] = circ_mean(eec_go{i}(hinds,:));
%     
    %PLOT Target
    subplot(1,length(DAYS),i); hold on;
    patch([xsforplot_target fliplr(xsforplot_target)],[dll_l.t fliplr(dul_l.t)],'b','FaceAlpha',0.5,'EdgeAlpha',0.5);
    patch([xsforplot_target fliplr(xsforplot_target)],[dll_h.t fliplr(dul_h.t)],'r','FaceAlpha',0.5,'EdgeAlpha',0.5);
    plot(xsforplot_target,dml.t,'b');
    plot(xsforplot_target,dmh.t,'r');
    
    %PLOT Go
    patch([xsforplot_go fliplr(xsforplot_go)],[dll_l.g fliplr(dul_l.g)],'b','FaceAlpha',0.5,'EdgeAlpha',0.5);
    patch([xsforplot_go fliplr(xsforplot_go)],[dll_h.g fliplr(dul_h.g)],'r','FaceAlpha',0.5,'EdgeAlpha',0.5);
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

%%
counts = DAYS_pred{1};
for i = 1:length(counts)
    
    spikecounts = counts{i};
    
    lcount(i,:) = sum(spikecounts(linds,:));
    hcount(i,:) = sum(spikecounts(hinds,:));
    
end

l_activ = sum(lcount);
h_activ = sum(hcount);

l_activ_boot = bootstrp(10000,@sum,lcount);
h_activ_boot = bootstrp(10000,@sum,hcount);

l_ord = sort(l_activ_boot);
h_ord = sort(h_activ_boot);

for i = 1:size(l_ord,2)
l_bottom(i) = l_ord(round(length(l_activ_boot(:,i))*0.025),i);
l_top(i) = l_ord(length(l_activ_boot(:,i)) - round(length(l_activ_boot(:,i))*0.025) + 1,i);

h_bottom(i) = h_ord(round(size(h_activ_boot(:,i),1)*0.025),i);
h_top(i) = h_ord(size(h_activ_boot(:,i),1) - round(size(h_activ_boot(:,i),1)*0.025) + 1,i);

end
times = (time_bins(2:end)+time_bins(1:end-1))./2;

figure; hold on;
patch([times fliplr(times)],[l_top fliplr(l_bottom)],'b','FaceAlpha',0.5,'EdgeAlpha',0.5);
patch([times fliplr(times)],[h_top fliplr(h_bottom)],'r','FaceAlpha',0.5,'EdgeAlpha',0.5);


