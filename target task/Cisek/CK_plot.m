brain_area = 'M1';

%%%% SPATIAL PRECISION %%%%%%%%%%%%%%%%%%%
spatial_cents = -pi:pi/16:pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% FOR PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%
fire_align = {6,7,8,9};%,7,9};
fire_ranges = {-200:50:1500, 0:50:500, 0:50:1500, 0:50:500};
BPDS = PDS;
tunes_used = {1,1,3,4};%,1,1,1};%,1,2};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FR_thresh = 0;
plotday = 2;
CO_index = 1;
reachdir_col = 17;
col_sort = 11;
make_symmetric = false;
do_peak_fit = false;

if exist('maxspeeds_co','var'); speed_input = maxspeeds_co; else speed_input = []; end
co_units = eval(sprintf('alldays(1).%s_units',brain_area));
if strcmp(brain_area,'M1'); co_units(end-1:end) = []; end

%% Tuning
tunereps = cell(length(fire_align),1);
for i = 1:length(fire_align)
    tunereps{i} = repmat(tunes_used{i},1,length(fire_ranges{i})-1);
end
tune_use = horzcat(tunereps{:});
%% Get Firing Rates
pred_days = [CO_index plotday];
[unit_counts,firing_absolute,spike_count_raw] = deal(cell(length(pred_days),1));
for day=1:length(pred_days)
    prediction_day = pred_days(day);
    for loopthrough = 1:length(fire_align)
        time_bins = fire_ranges{loopthrough};
        time_align = fire_align{loopthrough};

        day_pred = prediction_day;
        tt_pred = alldays(day_pred).tt;

        % find rasters for each unit
        unit_pred = cell(length(co_units),1);
        for q = 1:length(co_units)

            clc;
            fprintf('Trial block: %d/%d\nAlignment: %d/%d\nUnit: %d/%d\n',...
                day,length(pred_days),loopthrough,length(fire_align),q,length(co_units));

            [rast_out,rast_inds] = raster_plot(co_units{q},tt_pred,[time_bins(1)./1000 time_bins(end)./1000],...
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
            unit_counts{loopthrough}{q} = rast./repmat(diff(time_bins)/1000,size(rast,1),1);
        end
    end
    % 
    for loopthrough = 1:length(fire_align)
        time_bins = fire_ranges{loopthrough};
        time_align = fire_align{loopthrough};

        [trial_counts] = deal(cell(length(prediction_day),1));
        for bin = 1:length(time_bins) - 1
            clc; fprintf('Day: %d\nloop: %d/2\nbin %d/%d\n',day,loopthrough,bin,length(time_bins)-1);

            day_pred = prediction_day;
            tt_pred = alldays(day_pred).tt;  

            for i = 1:size(tt_pred,1)   
                for q = 1:length(co_units)       
                    trial_counts{bin}(i,q) = unit_counts{loopthrough}{q}(i,bin);
                end
            end
            firing_absolute{day}{loopthrough}{bin} = trial_counts{bin};
            spike_count_raw{day}{loopthrough}{bin} = trial_counts{bin}*(diff(time_bins(1:2))/1000);
        end

    end
    
end
FIRING = horzcat(firing_absolute{2}{:});
FIRING_CO = horzcat(firing_absolute{1}{:});
tt_fire = alldays(plotday).tt;
cofire = vertcat(FIRING_CO{:});
%% Likelihoods
likes = flipud(unique(tt_fire(:,col_sort)));
like_ind = cell(1,1);
for i = 1:length(likes)
    like_ind{1}{i} = find(tt_fire(:,col_sort)==likes(i));
end
%% Baselines
tbt_baselines = zeros(size(alldays(plotday).tt,1),length(co_units));
tpre2 = 0;
tpre1 = -200;
for i = 1:length(co_units) % Loop through units
    clc; fprintf('Calculating Baselines: (%d/%d)\n',i,length(co_units));

    [tbt_rast,tbt_inds] = raster_plot(co_units{i},alldays(plotday).tt,[tpre1 tpre2]./1000,6,0,'none');
    m = nansum(vertcat(tbt_rast{:}),2)./(size(tbt_rast{1},2)/1000);
%     tbt_baselines(vertcat(tbt_inds{:}),i) = m;
    tbt_baselines(vertcat(tbt_inds{:}),i) = nanmean(m);
end
%% DO DIRECTIONAL ORDERING
[DPD,GOODS] = deal(cell(length(FIRING),1));
for bin = 1:length(DPD)
    best_PDS = BPDS(:,tune_use(bin));
    repPD = repmat(best_PDS',size(tt_fire,1),1);
    reprch = repmat(tt_fire(:,reachdir_col),1,length(best_PDS));
    dfrompd = circ_dist(repPD,reprch);
    goods = ~isnan(best_PDS);
    DPD{bin} = dfrompd;
%     DPD{bin} = mod(dfrompd+pi/2,2*pi)-pi;%;{1}{end};
    GOODS{bin} = goods;
end
%dPD = mod(dfrompd(:,goods)+pi/2,2*pi)-pi;%;{1}{end};
space_size = diff(spatial_cents(1:2));
%% DO AVERAGING ACROSS TRIALS
[fire_trial,fire_all,dirs_all,Fire_dir,FIRE_dir,DS,FRS,BOOTDIR,BOOTALL,dpd_all,FRSt] = ...
    deal(cell(length(likes),1));
MEANBIN = zeros(length(likes),length(FIRING));
for lik = 1:length(likes)
    
    inds = find(tt_fire(:,col_sort)==likes(lik));
    for bin = 1:length(FIRING)
        goodmat = repmat(GOODS{bin}',length(inds),1);
 
        FIRE_norm_good = FIRING{bin}(inds,:).*goodmat - tbt_baselines(inds,:).*goodmat;
        DPD_good = DPD{bin}(inds,:).*goodmat;
        fire_all{lik}(:,bin) = reshape(FIRE_norm_good,[],1);
        dpd_all{lik}(:,bin) = reshape(DPD_good,[],1);
        fire_trial{lik}{bin} = FIRE_norm_good;

        clc; fprintf('averaging...\ncondition: %d/%d\ntime bin: %d/%d\n',lik,length(likes),bin,length(FIRING));

        for dir = 1:length(spatial_cents)
            if make_symmetric
                dirinds = find(abs(circ_dist(dpd_all{lik}(:,bin),abs(spatial_cents(dir))))<(space_size/2));
            else
                dirinds = find(abs(circ_dist(dpd_all{lik}(:,bin),spatial_cents(dir)))<(space_size/2));
            end
            Fire_dir{lik}(dir,bin) = nanmean(fire_all{lik}(dirinds,bin),1);
%             Fire_dir{lik}(dir,bin) = nanvar(fire_all{lik}(dirinds,bin),1);
%             FIRE_dir{lik}{dir}(:,bin) = fire_all{lik}(dirinds,bin);
        end 

        for trls = 1:size(fire_trial{lik}{bin},1)
            direcs = DPD{bin}(inds(trls),:);
            fires = fire_trial{lik}{bin}(trls,:);
            
            for dir = 1:length(spatial_cents)

                dirinds = find(abs(circ_dist(direcs,spatial_cents(dir)))<(space_size/2));
                
                FRSt{lik}{bin}(trls,dir) = nanmean(fires(dirinds));
            end

%                 DS{lik}(trls,:) = direcs;
%                 FRS{lik}(trls,bin,:) = fires;
        end
%         clc; fprintf('averaging...\ncondition: %d/%d\ntime bin: %d/%d\n',lik,length(likes),bin,length(FIRING)); 
    end
end
%% DO PLOTTING
breakpoint = fire_ranges{1}(end);
mincol = min(min(horzcat(Fire_dir{:})));
maxcol = max(max(horzcat(Fire_dir{:})));

time_cell = cellfun(@(x) (x(1:end-1)+x(2:end))/2,fire_ranges,'UniformOutput',0);
if length(time_cell)>1; 
    for i = 2:length(time_cell)
        time_cell{i} = time_cell{i}+max(time_cell{i-1}); 
    end
end
time_cents = horzcat(time_cell{:});

figure; %hold on; 
MC = [mincol maxcol];
TS = cell(length(likes),1);
for i = 1:length(likes)
    
    TS{i} = [Fire_dir{i}(3,:); mean(Fire_dir{i}([2 4],:)); Fire_dir{i}(1,:)];
%     TS{i} = [Fire_dir{i}(2,:); Fire_dir{i}(1,:)];
    
    subplot(length(likes),1,i); 

    imagesc(time_cents,spatial_cents,Fire_dir{i},MC); %[mincol,maxcol]); 
    hold on; 
    plot([0 0],[-1 1]*(pi+pi/64),'w','LineWidth',3);
    plot(breakpoint*[1 1],[-1 1]*(pi+pi/32),'w','LineWidth',3);
    %patch([800 1000 1000 800],[-1 -1 1 1]*(pi+pi/32),'w','EdgeColor','w');
%     plot(time_cents,MEANBIN(i,:),'w','LineWidth',3);
    eventtimes = cellfun(@(x) max(x),time_cell); eventtimes(end) = []; eventtimes = [0 eventtimes];
    
    box off
    title(sprintf('%.1f',unique(alldays(plotday).tt(alldays(plotday).tt(:,col_sort)==likes(i),col_sort))),'FontSize',16);
end
%% DO Difference
% figure;
% imagesc(time_cents,spatial_cents,Fire_dir{2}-Fire_dir{1});%,[mincol,maxcol]);
% hold on;
% plot([0 0],[-1 1]*(pi+pi/64),'w','LineWidth',3);
% plot(breakpoint*[1 1],[-1 1]*(pi+pi/32),'w','LineWidth',3);
% box off
% title(sprintf('%d',likes(i)),'FontSize',16);
%% DO Peak fitting
if do_peak_fit
    [decode_dir,boot_decode,decode_bounds] = deal(cell(length(FRS),1));
    figure; hold on; 
    clrs = {'b','r','c','m'};
    pvs = zeros(2,length(time_cents));
    peakform = zeros(2,1); 
    for i = 1:length(FRS) % Likelihood conditions

        for j = 1:size(FRS{i},2) % Time bins

            for k = 1:size(FRS{i},1) % Trials

                tbt_act = reshape(FRS{i}(k,j,:),[],1);
                tbt_dir = DS{i}(k,:)';

    %             b = glmfit([cos(tbt_dir) sin(tbt_dir)],tbt_act);
    %             
    %             decode_dir{i}(k,j) = abs(atan2(b(3),b(2)));
    %             

                repdirs = repmat(tbt_dir,1,length(spatial_cents));
                repbins = repmat(spatial_cents,size(tbt_dir,1),1);
                repacts = repmat(tbt_act,1,length(spatial_cents));
                dfspaces = circ_dist(repdirs,repbins);
                nearestbin = cell2mat(cellfun(@(x) abs(x)==min(abs(x)),num2cell(dfspaces,2),'UniformOutput',0));
                fin_bin = sum(repbins.*nearestbin,2);
                av_in_bin = sum(repacts.*nearestbin)./sum(nearestbin,1);

                decode_dir{i}(k,j) = circ_mean(spatial_cents(av_in_bin==max(av_in_bin))');

            end
            pvs(i,j) = circ_vtest(decode_dir{i}(:,j),0);
            clc; fprintf('%d - %d\n',i,j);

            [decode_bounds{i}(1,j),decode_bounds{i}(2,j),boot_decode{i}(:,j)] = ...
                boot_bounds(1000,@circ_mean,decode_dir{i}(:,j),2.5,97.5);
        end
        peakform(i) = find(pvs(i,:)<0.05,1,'first');
        plot(time_cents,circ_mean(boot_decode{i}),clrs{i}); 
        patch([time_cents fliplr(time_cents)],[decode_bounds{i}(1,:) fliplr(decode_bounds{i}(2,:))],...
            clrs{i},'FaceAlpha',0.25,'EdgeAlpha',0);
        plot(time_cents(peakform(i)).*[1 1],[-pi,pi],clrs{i},'LineWidth',3);
    end
end
 