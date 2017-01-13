brain_area = 'PMd';

loop_alignments = {'target','go'};
% loop_alignments = {'target'};
loop_ranges = {[-200:100:900],[-100:100:500]};
% loop_ranges = {[50 250]};

prebinnum = 2;

boot_type = 'trials';
CO_index = 1;
prediction_day = 6;
reachdir_col = 10;

%colors_of_plot = {'k','b','r','c'};
colors_of_plot = mat2cell(repmat(linspace(0,1,8)',1,3),ones(1,8),3);
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
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        baseline_reps = repmat(baseline_levels',length(tt_pred),1);
        
        from_baseline = trial_counts{bin}-baseline_reps;
        FR_changePD = nanmean(from_baseline.*PDinds,2);
        FR_changeOD = nanmean(from_baseline.*ODinds,2);
        FR_changeORTH = nanmean(from_baseline.*ORTHinds,2);
        FR_changeUT = nanmean(from_baseline.*UTinds,2);

        FR_change = FR_changeOD;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        set_of_inds = flipud(unique(alldays(prediction_day).tt(:,3)));
        like_ind = cell(length(set_of_inds),1);
        for qq = 1:length(set_of_inds)

            if t1 <= 0 && t2 <= 0 && strcmp(loop_alignments{loopthrough},'target')
                like_ind{qq} = find(ismember(tt_pred(:,3),set_of_inds));
            else
                like_ind{qq} = find(ismember(tt_pred(:,3),set_of_inds(qq)));
            end

            midbound(qq,bin) = mean(FR_change(like_ind{qq}));  
            mbound{loopthrough}(qq,bin) = midbound(qq,bin);

            if strcmp(boot_type,'trials')

                [lowerbound(qq,bin),upperbound(qq,bin)] = ...
                    boot_bounds(10000,@mean,FR_change(like_ind{qq}),2.5,97.5);     

            elseif strcmp(boot_type,'neurons')

                boot_across_neurs = mean(firing_diffs{bin}(like_ind{qq},:),1);

                [lowerbound(qq,bin),upperbound(qq,bin)] = ...
                    boot_bounds(10000,@mean,boot_across_neurs,2.5,97.5);    

            elseif strcmp(boot_type,'all')

                neuron_trials = firing_diffs{bin}(like_ind{qq},:);
                cat_neur_tri = reshape(neuron_trials,1,numel(neuron_trials));

                [lowerbound(qq,bin),upperbound(qq,bin)] = ...
                    boot_bounds(10000,@mean,cat_neur_tri,2.5,97.5);     
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