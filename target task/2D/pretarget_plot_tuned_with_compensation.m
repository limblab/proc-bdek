brain_area = 'PMd';

loop_alignments = {'target','go'};
% loop_alignments = {'target'};
loop_ranges = {[-200:100:900],[-100:100:500]};
% loop_ranges = {[50 250]};

prebinnum = 2;

boot_type = 'trials';
CO_index = 1;
prediction_day = 2;
reachdir_col = 10;

colors_of_plot = {'k','b','r','c'};
%colors_of_plot = mat2cell(repmat(linspace(0,1,8)',1,3),ones(1,8),3);
faceA = 1;
edgeA = 1;

% Initialize and precompute
dtime = unique(diff(loop_ranges{1}));
[mbound,totd] = deal(cell(1,length(loop_alignments)));
co_units = eval(sprintf('alldays(1).%s_units',brain_area));
if strcmp(brain_area,'M1'); co_units(end-1:end) =[]; end
tt_pred = alldays(prediction_day).tt;

%% Calculate Baseline Levels
baseraster = cell(length(co_units),1);
tpre2 = 0;
tpre1 = 0-prebinnum*dtime;
for i = 1:length(co_units) % Loop through units
    clc; fprintf('Calculating Baselines: (%d/%d)\n',i,length(co_units));

    baseraster{i} = raster_plot(co_units{i},alldays(CO_index).tt,[tpre1 tpre2]./1000,'target',0,'none');
end
baseline_levels = cell2mat(cellfun(@(x) nanmean(nansum(x{1},2)./(size(x{1},2)/1000),1),baseraster,'UniformOutput',0));   

%% Do PD calculations

reppd = repmat(best_PDS',size(alldays(prediction_day).tt,1),1);
reprch = repmat(alldays(prediction_day).tt(:,reachdir_col),1,size(reppd,2));

dfrompd = circ_dist(reppd,reprch);

PDinds = double(abs(circ_dist(reppd,reprch)) <= pi/4); PDinds(PDinds==0)=NaN;
ODinds  = double(abs(circ_dist(reppd,reprch)) >= 3*pi/4); ODinds(ODinds==0)=NaN;
ORTHinds = double(abs(circ_dist(reppd,reprch))>pi/4 & abs(circ_dist(reppd,reprch))<3*pi/4); ORTHinds(ORTHinds==0)=NaN;
UTinds = double(isnan(reppd)); UTinds(UTinds==0)=NaN;

%OPranges = repmat(op_ranges',size(reprch,1),1); OPranges(OPranges==0) = NaN;

%% Do bin-by-bin tuning for reach direction compensation
[act_predict,expect_counts] = deal(cell(length(loopthrough),1));
for loopthrough = 1:length(loop_alignments)  
tune_align = loop_alignments{loopthrough};
tune_ranges = loop_ranges{loopthrough};
    for ranges = 1:length(tune_ranges)-1
        tune_range = tune_ranges(ranges:ranges+1);
        for i = 1:length(co_units) % Loop through units
            clc; fprintf('%d/%d Tuning: (%d/%d)\n',ranges,length(tune_ranges)-1,i,length(co_units));

            % set time bin edges
            t1 = tune_range(1);
            t2 = tune_range(2);

            if  t1 <= 0 && t2 <= 0 && loopthrough == 1
                [cents,tune] = co_tuning(co_units,[nan; alldays(CO_index).tt(1:end-1,reachdir_col)],...
                alldays(CO_index).tt,t1,t2,tune_align,i);                
            else 
                [cents,tune] = co_tuning(co_units,alldays(CO_index).tt(:,10),...
                alldays(CO_index).tt,t1,t2,tune_align,i);
            end

            tune_cents = [cents-2*pi cents cents+2*pi];
            tune_pad = [tune;tune;tune];
            
            interp_cents = tune_cents(1):.01:tune_cents(end);
            tune_interp = interp1(tune_cents,tune_pad,interp_cents);
            
            wrapped_cents = interp_cents(interp_cents >=0 & interp_cents< 2*pi);
            wrapped_tune = tune_interp(interp_cents >= 0 & interp_cents <2*pi); 
           
            act_predict{loopthrough}{ranges}(i,:) = wrapped_tune;
        end
    end 
end
move_index = zeros(size(tt_pred,1),1);
for i = 1:size(tt_pred,1)
    dind = abs(circ_dist(tt_pred(i,10),wrapped_cents));
    move_index(i) = find(dind==min(dind));
end
%% Get Counts
unit_counts = cell(length(loop_alignments),1);
for loopthrough = 1:length(loop_alignments)
    
    tune_align = loop_alignments{loopthrough};
    tune_ranges = loop_ranges{loopthrough};
    time_bins = tune_ranges;
    time_align = tune_align;

    % find rasters for each unit
    [unit_pred] = deal(cell(length(co_units),1));
    for q = 1:length(co_units)

        clc; fprintf('Computing Rasters...\nAlignment: %d/%d\nUnit: %d/%d\n',...
            loopthrough,length(loop_alignments),q,length(co_units));

        [rast_out,rast_inds] = raster_plot(co_units{q},tt_pred,[time_bins(1)./1000 time_bins(end)./1000],...
            time_align,0,'none');

        out = nan(size(vertcat(rast_out{:})));
        for m = 1:length(rast_out)
            out(rast_inds{m},:) = rast_out{m};
        end
        out(isnan(out))=0;

        rast = zeros(size(out,1),length(time_bins)-1);
        binsizes = time_bins - time_bins(1); binsizes(1) = 0;
        for v = 1:length(time_bins)-1
            rast(:,v) = sum(out(:,(binsizes(v)+1):binsizes(v+1)),2);
        end
        unit_counts{loopthrough}{q} = rast./repmat(diff(time_bins)/1000,size(rast,1),1); % in FIRING RATE
    end
    
end

%% Calculate Metrics and Plot
figure; hold on;
for loopthrough = 1:length(loop_alignments)
    
    tune_align = loop_alignments{loopthrough};
    tune_ranges = loop_ranges{loopthrough};
    time_bins = tune_ranges;
    time_align = tune_align;
    
    [lowerbound, upperbound,midbound] = ...
        deal(zeros(length(unique(alldays(prediction_day).tt(:,3))),length(time_bins)-1));
    trial_counts = cell(length(time_bins)-1,1);
    
    for bin = 1:length(time_bins) - 1

        tune_range = tune_ranges(bin:bin+1);
        t1 = tune_range(1);
        t2 = tune_range(2);

        clc; fprintf('Calculating for each time bin...\nAlignment: %d/%d\nBin %d/%d\n',...
            loopthrough,length(loop_alignments),bin,length(time_bins)-1);

        for i = 1:size(tt_pred,1)   
            for q = 1:length(co_units)       
                trial_counts{bin}(i,q) = unit_counts{loopthrough}{q}(i,bin); % in FIRING RATE
            end
            expect_counts{loopthrough}{bin}(i,:) = act_predict{loopthrough}{bin}(:,move_index(i))';
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        baseline_reps = repmat(baseline_levels',length(tt_pred),1);
        
        from_baseline = trial_counts{bin}-baseline_reps;
        from_zerocond = trial_counts{bin}-expect_counts{loopthrough}{bin};
        
        %%%
        voi = from_zerocond;
        %%%
        
        FR_changePD = voi.*PDinds;
        FR_changeOD = voi.*ODinds;
        FR_changeORTH = voi.*ORTHinds;
        FR_changeUT = voi.*UTinds;
        
        %%%%%%%%%%%%%%%%%%%%%%\
        voif = FR_changeOD; %%|
        %%%%%%%%%%%%%%%%%%%%%%/
        
        FR_change = nanmean(voif,2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set_of_inds = flipud(unique(alldays(prediction_day).tt(:,3)));
        like_ind = cell(length(set_of_inds),1);
        for qq = 1:length(set_of_inds)

            if t1 <= 0 && t2 <= 0 && strcmp(loop_alignments{loopthrough},'target')
                like_ind{qq} = find(ismember(tt_pred(:,3),set_of_inds));
            else
                like_ind{qq} = find(ismember(tt_pred(:,3),set_of_inds(qq)));
            end
            
%             midbound(qq,bin) = nanmean(FR_change(like_ind{qq}));  
            if strcmp(boot_type,'trials')

                midbound(qq,bin) = nanmean(reshape(voif(like_ind{qq},:),[],1));
                mbound{loopthrough}(qq,bin) = midbound(qq,bin);
                
                [lowerbound(qq,bin),upperbound(qq,bin)] = ...
                    boot_bounds(10000,@mean,FR_change(like_ind{qq}),2.5,97.5);     

            elseif strcmp(boot_type,'neurons')

                boot_across_neurs = nanmean(voif(like_ind{qq},:),1);
                boot_across_nonnan = boot_across_neurs(~isnan(boot_across_neurs));

                midbound(qq,bin) = nanmean(boot_across_neurs);
                mbound{loopthrough}(qq,bin) = midbound(qq,bin);
                
                [lowerbound(qq,bin),upperbound(qq,bin)] = ...
                    boot_bounds(10000,@nanmean,boot_across_nonnan,2.5,97.5);    

            elseif strcmp(boot_type,'all')
                
                boot_across_neurs = nanmean(voif(like_ind{qq},:),1);
                boot_across_nonnan = boot_across_neurs(~isnan(boot_across_neurs));
                midbound(qq,bin) = nanmean(boot_across_neurs);
                
                sem_all = std(reshape(voif(~isnan(voif)),[],1))./sqrt((sum(sum(~isnan(voif)))));

                lowerbound(qq,bin) = midbound(qq,bin) - sem_all;
                upperbound(qq,bin) = midbound(qq,bin) + sem_all;
  
            else

                fprintf('Improper bootstrap dimension: Bootstrapping across all\n'); 

                neuron_trials = firing_diffs{bin}(like_ind{qq},:);
                cat_neur_tri = reshape(neuron_trials,1,numel(neuron_trials));

                [lowerbound(qq,bin),upperbound(qq,bin)] = ...
                    boot_bounds(10000,@mean,cat_neur_tri,2.5,97.5);   

            end

            totd{loopthrough}{qq}(:,bin) = [midbound(qq,bin); lowerbound(qq,bin); upperbound(qq,bin)];
        end
    end
    %
    dayind = 1;

    if loopthrough == 1

        xsforplotT = .5*(time_bins(1:end-1)+time_bins(2:end));
        max_targ_x = max(xsforplotT);
        for i = 1:length(unique(alldays(prediction_day).tt(:,3)))
            plot(xsforplotT,midbound(i,:),'Color',colors_of_plot{i});
            patch([xsforplotT fliplr(xsforplotT)],[lowerbound(i,:) fliplr(upperbound(i,:))],...
                colors_of_plot{i},'FaceAlpha',faceA,'EdgeAlpha',edgeA);
        end
        xsfp{1} = xsforplotT;

    else

        xsatgo = .5*(time_bins(1:end-1)+time_bins(2:end));
        xsforplotG = xsatgo + max_targ_x + 100;
        min_go_x = min(xsatgo);

        for i = 1:length(unique(alldays(prediction_day).tt(:,3)))
            plot(xsforplotG,midbound(i,:),'Color',colors_of_plot{i});
            patch([xsforplotG fliplr(xsforplotG)],[lowerbound(i,:) fliplr(upperbound(i,:))],...
                colors_of_plot{i},'FaceAlpha',faceA,'EdgeAlpha',edgeA);
        end 
        xsfp{2} = xsforplotG;
    end
end
%
FDc = horzcat(FD{:});
[pos_totd] = deal(cell(1,2));
for prepost = 1:length(loop_alignments)
    positiveinds = xsfp{prepost} > 0;
    pos_totd{prepost} = cellfun(@(x) x(:,positiveinds), totd{prepost},'UniformOutput',0);
end

plot(max_targ_x*[1 1],[-.5 1.5],'k--');
plot([0 0],[-.5 1.5],'k--');
plot((max_targ_x+min_go_x+100)*[1 1],[-.5 1.5],'k--');
xlabel('Time From Target On (ms)','FontSize',14);
ylabel('Change in FR (% modulation)','FontSize',14);
title(sprintf('%s (%s) - %s/%s/%s',monkey,brain_area,MO,DA,YE),'FontSize',16);