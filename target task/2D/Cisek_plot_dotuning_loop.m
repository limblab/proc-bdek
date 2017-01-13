brain_area = BRAIN_AREA;

%%%% SPATIAL PRECISION %%%%%%%%%%%%%%%%%%%
spatial_cents = -pi:pi/8:pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% FOR CALCULATING TUNING CURVES %%%%%%%%
% loop_alignments = {'target'};%,'go'};
% loop_ranges = {[-400 0]};%,[50 250]};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FR_thresh = 0;
plotday = 2;
CO_index = 1;
reachdir_col = 10;
col_sort = 3;
make_symmetric = false;
do_peak_fit = true;

if exist('maxspeeds_co','var'); speed_input = maxspeeds_co; else speed_input = []; end
co_units = eval(sprintf('alldays(1).%s_units',brain_area));
%% DIRECTIONAL TUNING
if ~exist('best_PDS','var') 
    best_PDS = tuning_types_PD_func(alldays,brain_area,CO_index,reachdir_col,loop_alignments,loop_ranges);
end
%% Get Firing Rates
pred_days = [CO_index plotday];
[unit_counts,firing_absolute,spike_count_raw] = deal(cell(length(pred_days),1));
for day=1:length(pred_days)
    prediction_day = pred_days(day);
    for loopthrough = 1:length(fire_align)
        time_bins = fire_ranges{loopthrough};
        time_align = fire_align{loopthrough};

        day_units = eval(sprintf('alldays(1).%s_units',brain_area));
        if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end

        day_pred = prediction_day;
        tt_pred = alldays(day_pred).tt;

        % find rasters for each unit
        unit_pred = cell(length(day_units),1);
        for q = 1:length(day_units)

            clc;
            fprintf('Trial block: %d/%d\nAlignment: %d/%d\nUnit: %d/%d\n',...
                day,length(pred_days),loopthrough,length(fire_align),q,length(day_units));

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

            day_units = eval(sprintf('alldays(1).%s_units',brain_area));
            if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end

            for i = 1:size(tt_pred,1)   
                for q = 1:length(day_units)       
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
if do_baselines
    tbt_baselines = zeros(size(alldays(plotday).tt,1),length(co_units));
    tpre2 = 0;
    tpre1 = -200;
    for i = 1:length(co_units) % Loop through units
        clc; fprintf('Calculating Baselines: (%d/%d)\n',i,length(co_units));

        [tbt_rast,tbt_inds] = raster_plot(co_units{i},alldays(plotday).tt,[tpre1 tpre2]./1000,'target',0,'none');
        m = nansum(vertcat(tbt_rast{:}),2)./(size(tbt_rast{1},2)/1000);
    %     tbt_baselines(vertcat(tbt_inds{:}),i) = m;
        tbt_baselines(vertcat(tbt_inds{:}),i) = nanmean(m);
    end
end
%% DO DIRECTIONAL ORDERING
repPD = repmat(best_PDS',size(tt_fire,1),1);
reprch = repmat(tt_fire(:,reachdir_col),1,length(best_PDS));

dfrompd = circ_dist(repPD,reprch);

goods = ~isnan(best_PDS);
dPD = dfrompd(:,goods);
%dPD = mod(dfrompd(:,goods)+pi/2,2*pi)-pi;%;{1}{end};
space_size = diff(spatial_cents(1:2));
%% DO AVERAGING ACROSS TRIALS
[fire_trial,fire_all,dirs_all,Fire_dir,FIRE_dir,DS,FRS,BOOTDIR,BOOTALL] = ...
    deal(cell(length(likes),1));
MEANBIN = zeros(length(likes),length(FIRING));
for lik = 1:length(likes)
    
    inds = find(tt_fire(:,col_sort)==likes(lik));
    dirs_all{lik} = reshape(dPD(inds,:),[],1);
    for bin = 1:length(FIRING)
 
        FIRE_norm_good = FIRING{bin}(inds,goods) - tbt_baselines(inds,goods);
        fire_all{lik}(:,bin) = reshape(FIRE_norm_good,[],1);
        fire_trial{lik}{bin} = FIRE_norm_good;

        clc; fprintf('averaging...\ncondition: %d/%d\ntime bin: %d/%d\n',lik,length(likes),bin,length(FIRING));

%         MEANBIN(lik,bin) = glmhelper([dirs_all{lik} fire_all{lik}(:,bin)]);
%         MEANBIN(lik,bin) = mean(fire_all{lik}(dirs_all{lik}>0,bin)) / ...
%                            mean(fire_all{lik}(dirs_all{lik}<0,bin));
    end
    
    for dir = 1:length(spatial_cents)
        if make_symmetric
            dirinds = find(abs(circ_dist(dirs_all{lik},abs(spatial_cents(dir))))<(space_size/2));
        else
            dirinds = find(abs(circ_dist(dirs_all{lik},spatial_cents(dir)))<(space_size/2));
        end
        Fire_dir{lik}(dir,:) = nanmean(fire_all{lik}(dirinds,:),1);
        FIRE_dir{lik}{dir} = fire_all{lik}(dirinds,:);
    end 

    for bin = 1:length(FIRING)
        for trls = 1:size(fire_trial{lik}{bin},1)
            
            direcs = dPD(inds(trls),:);
            fires = fire_trial{lik}{bin}(trls,:);
            
%             [BOOTDIR{lik}{trls}(1,bin),BOOTDIR{lik}{trls}(2,bin),~,BOOTALL{lik}{trls}(:,bin)] = ...
%             boot_bounds(100,@glmhelper,[direcs',fires'],2.5,97.5);
           
%             for bstp = 1:1000
%                 [bb] = glmfit([cos(direcs') sin(direcs')],fires);
%                 BOOTDIR{lik}(trls,bstp) = atan2(bb(3),bb(2));
%                 clc; fprintf('%d %d %d',bin,trls,bstp);
%             end
            DS{lik}(trls,:) = direcs;
            FRS{lik}(trls,bin,:) = fires;
        end
        clc; fprintf('averaging...\ncondition: %d/%d\ntime bin: %d/%d\n',lik,length(likes),bin,length(FIRING)); 
    end
end
%% DO PLOTTING
breakpoint = fire_ranges{1}(end);
mincol = min(min(horzcat(Fire_dir{:})));
maxcol = max(max(horzcat(Fire_dir{:})));

if ~exist('MC','var'); MC = [mincol maxcol]; end

time_cell = cellfun(@(x) (x(1:end-1)+x(2:end))/2,fire_ranges,'UniformOutput',0);
if length(time_cell)>1; time_cell{2} = time_cell{2}+breakpoint; end
time_cents = horzcat(time_cell{:});

% figure; hold on; 
% for i = 1:length(likes)
%     
%     subplot(length(likes),1,i); 
% %     figure; hold on; 
%     imagesc(time_cents,spatial_cents,Fire_dir{i},MC); %[mincol,maxcol]);
%     
%     hold on;
%     plot([0 0],[-1 1]*(pi+pi/64),'w','LineWidth',3);
%     plot(breakpoint*[1 1],[-1 1]*(pi+pi/32),'w','LineWidth',3);
%     %patch([800 1000 1000 800],[-1 -1 1 1]*(pi+pi/32),'w','EdgeColor','w');
% %     plot(time_cents,MEANBIN(i,:),'w','LineWidth',3);
%     
%     box off
%     title(sprintf('%d',likes(i)),'FontSize',16);
% end
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
    [decode_dir,boot_decode,decode_bounds,avs_in_bins,decode_unc,decode_std,decodeVAF] = deal(cell(length(FRS),1));
%     figure; hold on; 
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
%                 fin_bin = sum(repbins.*nearestbin,2);
                av_in_bin = sum(repacts.*nearestbin,1)./sum(nearestbin,1);
                avs_in_bins{i}{j}(k,:) = av_in_bin;
                
                nonempt = find(~isnan(av_in_bin));
                A = [ones(length(nonempt),1) cos(spatial_cents(nonempt)') sin(spatial_cents(nonempt)')];
                b = A\av_in_bin(nonempt)'; 
                decode_dir{i}(k,j) = atan2(b(3),b(2));
                decode_unc{i}(k,j) = sqrt(b(2).^2+b(3).^2);
                
%                 decode_dir{i}(k,j) = circ_mean(spatial_cents(av_in_bin==max(av_in_bin))');

            end
            pvs(i,j) = circ_vtest(decode_dir{i}(:,j),0);
            clc; fprintf('%d - %d\n',i,j);

            [decode_bounds{i}(1,j),decode_bounds{i}(2,j),boot_decode{i}(:,j)] = ...
                boot_bounds(1000,@circ_mean,decode_dir{i}(:,j),2.5,97.5);
        end
%         peakform(i) = find(pvs(i,:)<0.05,1,'first');
%         plot(time_cents,circ_mean(boot_decode{i}),clrs{i}); 
%         patch([time_cents fliplr(time_cents)],[decode_bounds{i}(1,:) fliplr(decode_bounds{i}(2,:))],...
%             clrs{i},'FaceAlpha',0.25,'EdgeAlpha',0);
%         plot(time_cents(peakform(i)).*[1 1],[-pi,pi],clrs{i},'LineWidth',3);
        
%           tbt_dir_nonan = tbt_dir; tbt_dir_nonan(isnan(tbt_dir_nonan)) = [];
          decode_dir_nonan = decode_dir{i}; decode_dir_nonan(isnan(decode_dir_nonan)) = [];
          
          decodeVAF{i}(1,:) = (circ_var(tbt_dir)-circ_var(decode_dir_nonan))./circ_var(tbt_dir);
          decode_std{i}(1,:) = circ_std(decode_dir_nonan);
%         [decode_std{i}(2,:),decode_std{i}(3,:)] = boot_bounds(1000,@circ_std,decode_dir{i},2.5,97.5);
%         plot(time_cents,decode_std{i}(1,:),clrs{i});
%         patch([time_cents fliplr(time_cents)],[decode_std{i}(2,:) fliplr(decode_std{i}(3,:))],clrs{i},'FaceAlpha',0.25,'EdgeColor','none');
    end
end
 