brain_area = 'PMd';

loop_alignments = {'target'};
% loop_alignments = {12};
% loop_ranges = {[-200:50:850],[-50:50:450]};
% loop_ranges = {[-200:100:800],[-50:100:400]};
loop_ranges = {[300 700]};%,[-100:100:400]};
% loop_ranges = {[-225:25:725]};
% loop_ranges = {[500 700]};

loop_alignmentsPD = {'target','go'};
% loop_alignmentsPD = {'target'};
loop_rangesPD = {[500 700],[50 250]};
% loop_rangesPD = {[600 800]};

baseline_window = 200;
boot_type = 'all';
do_expected = true;
use_speeds = false;

CO_index = 1;
prediction_day = 2;
reachdir_col = 10;
boot_number = 1000;

colors_of_plot = {'b','r','c','m'};
%colors_of_plot = mat2cell(repmat(linspace(0,1,8)',1,3),ones(1,8),3);
faceA = .25;
edgeA = 1;

% Initialize and precompute
dtime = unique(diff(loop_ranges{1}));
[mbound,lbound,ubound,totd] = deal(cell(1,length(loop_alignments)));
co_units = eval(sprintf('alldays(1).%s_units',brain_area));
if strcmp(brain_area,'M1'); co_units(end-1:end) =[]; end
tt_pred = alldays(prediction_day).tt;
coangs = alldays(CO_index).tt(:,reachdir_col);
sessangs = alldays(prediction_day).tt(:,reachdir_col);
%% Get preferred directions
if ~exist('best_PDS','var')
    best_PDS = tuning_types_PD_func(alldays,brain_area,CO_index,reachdir_col,loop_alignmentsPD,loop_rangesPD);
end
%% Calculate Average Baseline Levels
baseraster = cell(length(co_units),1);
for i = 1:length(co_units) % Loop through units
    clc; fprintf('Calculating Baselines: (%d/%d)\n',i,length(co_units));

    baseraster{i} = raster_plot(co_units{i},alldays(CO_index).tt,[-baseline_window,0]./1000,'target',0,'none');
end
baseline_levels = cell2mat(cellfun(@(x) nanmean(nansum(x{1},2)./(size(x{1},2)/1000),1),baseraster,'UniformOutput',0));  
tbt_baselines = repmat(baseline_levels',size(alldays(prediction_day).tt,1),1);
%% Calculate Trial-by-trial Baseline Levels
% tbt_baselines = zeros(size(alldays(prediction_day).tt,1),length(co_units));
% for i = 1:length(co_units) % Loop through units
%     clc; fprintf('Calculating Baselines: (%d/%d)\n',i,length(co_units));
%     
%     rast = raster_get(co_units{i},alldays(prediction_day).tt,[-baseline_window,0]./1000,'target');
%     tbt_baselines(:,i) = nansum(rast,2)./(baseline_window/1000);
% end
%% Do PD calculations

reppd = repmat(best_PDS',size(alldays(prediction_day).tt,1),1);
reprch = repmat(alldays(prediction_day).tt(:,reachdir_col),1,size(reppd,2));

dfrompd = circ_dist(reppd,reprch);

PDinds = double(abs(circ_dist(reppd,reprch)) <= pi/4); PDinds(PDinds==0)=NaN;
% PDinds = double(abs(circ_dist(reppd,reprch)) <= pi/2); PDinds(PDinds==0)=NaN;
ODinds  = double(abs(circ_dist(reppd,reprch)) >= 3*pi/4); ODinds(ODinds==0)=NaN;
% ODinds  = double(abs(circ_dist(reppd,reprch)) >= pi/2); ODinds(ODinds==0)=NaN;
ORTHinds = double(abs(circ_dist(reppd,reprch))>pi/4 & abs(circ_dist(reppd,reprch))<3*pi/4); ORTHinds(ORTHinds==0)=NaN;
UTinds = double(isnan(reppd)); UTinds(UTinds==0)=NaN;
ALLinds = double(~isnan(reppd)); ALLinds(ALLinds==0)=NaN;

%OPranges = repmat(op_ranges',size(reprch,1),1); OPranges(OPranges==0) = NaN;
%% Do bin-by-bin tuning for reach direction compensation
if do_expected
    
    if use_speeds
        speeds_co = kin_exam(BDF,alldays(CO_index).tt,[],[],[6 7]);
        maxspeeds_co = cellfun(@(x) max(x),speeds_co);
        speeds_pred = kin_exam(BDF,tt_pred,[],[],[6 7]);
        maxspeeds_pred = cellfun(@(x) max(x),speeds_pred);
    end

    [expect_counts,TC] = deal(cell(length(loop_alignments),1));
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
                
                [cents,tune] = co_tuning(co_units,alldays(CO_index).tt(:,reachdir_col),...
                                         alldays(CO_index).tt,t1,t2,tune_align,i);
                cents_rep = [cents-2*pi cents cents+2*pi]; tune_rep = [tune' tune' tune'];
                expect_counts{loopthrough}{ranges}(:,i) = interp1(cents_rep,tune_rep,sessangs);
                TC{loopthrough}{ranges}(:,i) = interp1(cents_rep,tune_rep,0:0.01:2*pi);
            end
        end 
    end
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

        out = raster_get(co_units{q},tt_pred,time_bins([1,end])./1000,time_align);
%         out(isnan(out))=0;

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
[mids,lows,highs,voif_all,rbs,av_all] = deal(cell(4,1));
tunelines = {'.-','s-','o-','--'};
for tuneloop = 1:4
    indexer = 0;
    randboot = cell(1,length(loop_alignments));
    for loopthrough = 1:length(loop_alignments)

        tune_align = loop_alignments{loopthrough};
        tune_ranges = loop_ranges{loopthrough};
        time_bins = tune_ranges;
        time_align = tune_align;

        [lowerbound, upperbound,midbound] = ...
            deal(zeros(length(unique(alldays(prediction_day).tt(:,3))),length(time_bins)-1));
        trial_counts = cell(length(time_bins)-1,1);

        for bin = 1:length(time_bins) - 1
            indexer = indexer+1;
            tune_range = tune_ranges(bin:bin+1);
            t1 = tune_range(1);
            t2 = tune_range(2);

            clc; fprintf('Calculating for each time bin...\nloop: %d\nAlignment: %d/%d\nBin %d/%d\n',...
                tuneloop,loopthrough,length(loop_alignments),bin,length(time_bins)-1);

            for i = 1:size(tt_pred,1)   
                for q = 1:length(co_units)       
                    trial_counts{bin}(i,q) = unit_counts{loopthrough}{q}(i,bin); % in FIRING RATE
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             baseline_reps = repmat(baseline_levels',length(tt_pred),1);
            
%             from_zerocond = trial_counts{bin}-expect_counts{loopthrough}{bin};
%             pred_baseline = (expect_counts{loopthrough}{bin}.*1000./(t2-t1))-baseline_reps;
%             expect_tbt_baseline = expect_counts{loopthrough}{bin}-tbt_baselines;
            from_tbt_baseline = trial_counts{bin}-tbt_baselines;
            %-------------------|
%             voi = expect_counts{loopthrough}{bin}-tbt_baselines;%.*1000./(t2-t1);%from_tbt_baseline;
%             voi = trial_counts{bin};
  
%             voi = from_zerocond;
            voi = from_tbt_baseline;
%             voi = trial_counts{bin}-(expect_counts{loopthrough}{bin});%.*1000./(t2-t1));
%             voi = expect_tbt_baseline;
            %-------------------|
            voi(voi>500) = NaN;

            FR_changePD = voi.*PDinds;
            FR_changeOD = voi.*ODinds;
            FR_changeORTH = voi.*ORTHinds;
            FR_changeUT = voi.*UTinds;
            FR_changeNT = voi.*(ODinds==1|UTinds==1);
            FR_changeTUNED = voi.*(PDinds==1|ODinds==1|ORTHinds==1);
            FR_changeALL = voi.*ALLinds;
            
            excitatory = double(voi>0); excitatory(~excitatory)=NaN;
            inhibitory = double(voi<0); inhibitory(~inhibitory)=NaN;
            
%             expect_depth = nanmean(voie.*PDinds,2)-nanmean(voie.*ODinds,2);
%             actual_depth = nanmean(voi.*PDinds,2)-nanmean(voi.*ODinds,2);
%             EX_IN = {excitatory,inhibitory};    
%             PD_OD = {FR_changePD, FR_changeOD};
            
            PD_OD = {FR_changePD, FR_changeOD, FR_changeORTH, FR_changeALL};%, voie.*PDinds, voie.*ODinds};
            AV_BIN = {trial_counts{bin}.*PDinds, trial_counts{bin}.*ODinds,trial_counts{bin}.*ORTHinds,trial_counts{bin}.*ALLinds};
            %%%%%%%%%%%%%%%%%%%%%%\
%             voif = actual_depth - expect_depth; %PD_OD{tuneloop}; %%|

            voif = PD_OD{tuneloop};
%             voif = nanvar(FR_changeALL,[],1);
%             voif = FR_changeALL;
%             voif = nanmean(PD_OD{1},2)-nanmean(PD_OD{2},2);

%             voif = nanvar(PD_OD{tuneloop},[],2)./nanmean(PD_OD{tuneloop},2);
%             voif = P{1};
            %%%%%%%%%%%%%%%%%%%%%%/
            av_all{tuneloop}{indexer} = AV_BIN{tuneloop};
            voif_all{tuneloop}{indexer} = voif;
            
            FR_change = nanmean(voif,2);
            FR_change(sum(isnan(voif),2)==size(voif,2)) = NaN;
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

                    %midbound(qq,bin) = nanmean(reshape(voif(like_ind{qq},:),[],1));
                    midbound(qq,bin) = nanmean(FR_change(like_ind{qq}));

                    [lowerbound(qq,bin),upperbound(qq,bin),~,randboot{loopthrough}{qq,bin}] = ...
                        boot_bounds(boot_number,@nanmean,FR_change(like_ind{qq}),2.5,97.5);     

                elseif strcmp(boot_type,'neurons')

                    boot_across_neurs = nanmean(voif(like_ind{qq},:),1);
                    boot_across_nonnan = boot_across_neurs(~isnan(boot_across_neurs));

                    midbound(qq,bin) = nanmean(boot_across_neurs);
                    mbound{loopthrough}(qq,bin) = midbound(qq,bin);


                    [lowerbound(qq,bin),upperbound(qq,bin)] = ...
                        boot_bounds(boot_number,@nanmean,boot_across_nonnan,2.5,97.5);    

                elseif strcmp(boot_type,'all')

                    bootcat = reshape(voif(like_ind{qq},:),[],1);
                    boot_nonan = bootcat(~isnan(bootcat));
                    midbound(qq,bin) = mean(boot_nonan);
                    
                    [lowerbound(qq,bin),upperbound(qq,bin)] = boot_bounds(1000,@mean,boot_nonan,2.5,97.5);
%                     
%                     boot_across_neurs = nanmean(voif(like_ind{qq},:),1);
%                     boot_across_nonnan = boot_across_neurs(~isnan(boot_across_neurs));
%                     midbound(qq,bin) = nanmean(boot_across_neurs);
% 
%                     sem_all = std(reshape(voif(~isnan(voif)),[],1))./sqrt((sum(sum(~isnan(voif)))));
% 
%                     lowerbound(qq,bin) = midbound(qq,bin) - sem_all;
%                     upperbound(qq,bin) = midbound(qq,bin) + sem_all;

                else

                    fprintf('Improper bootstrap dimension: Bootstrapping across all\n'); 

                    neuron_trials = firing_diffs{bin}(like_ind{qq},:);
                    cat_neur_tri = reshape(neuron_trials,1,numel(neuron_trials));

                    [lowerbound(qq,bin),upperbound(qq,bin)] = ...
                        boot_bounds(10000,@mean,cat_neur_tri,2.5,97.5);   

                end
                mbound{loopthrough}(qq,bin) = midbound(qq,bin);
                lbound{loopthrough}(qq,bin) = lowerbound(qq,bin);
                ubound{loopthrough}(qq,bin) = upperbound(qq,bin);

                totd{loopthrough}{qq}(:,bin) = [midbound(qq,bin); lowerbound(qq,bin); upperbound(qq,bin)];
            end
        end
        %
        dayind = 1;

        if loopthrough == 1

            xsforplotT = .5*(time_bins(1:end-1)+time_bins(2:end));
            max_targ_x = max(xsforplotT);
            for i = 1:length(unique(alldays(prediction_day).tt(:,3)))
                plot(xsforplotT,midbound(i,:),tunelines{tuneloop},'Color',colors_of_plot{i});

                patch([xsforplotT fliplr(xsforplotT)],[lowerbound(i,:) fliplr(upperbound(i,:))],...
                    colors_of_plot{i},'FaceAlpha',faceA,'EdgeAlpha',edgeA);

            end
            xsfp{1} = xsforplotT;

        else

            xsatgo = .5*(time_bins(1:end-1)+time_bins(2:end));
            xsforplotG = xsatgo + 1000;% max_targ_x + 100;
            min_go_x = min(xsatgo);

            for i = 1:length(unique(alldays(prediction_day).tt(:,3)))
                plot(xsforplotG,midbound(i,:),tunelines{tuneloop},'Color',colors_of_plot{i});

                patch([xsforplotG fliplr(xsforplotG)],[lowerbound(i,:) fliplr(upperbound(i,:))],...
                    colors_of_plot{i},'FaceAlpha',faceA,'EdgeAlpha',edgeA);

            end 
            xsfp{2} = xsforplotG;
        end
    end

    mids{tuneloop} = horzcat(mbound{:});
    lows{tuneloop} = horzcat(lbound{:});
    highs{tuneloop} = horzcat(ubound{:});
    rbs{tuneloop} = randboot; 
end

if exist('min_go_x','var')
    plot(max_targ_x*[1 1],[-.5 1.5],'k--');
    plot([0 0],[-.5 1.5],'k--');
    plot((max_targ_x+min_go_x+100)*[1 1],[-.5 1.5],'k--');
    xlabel('Time From Target On (ms)','FontSize',14);
    ylabel('Change in FR (% modulation)','FontSize',14);
    title(sprintf('%s (%s) - %s/%s/%s',monkey,brain_area,MO,DA,YE),'FontSize',16);
end
%%
% figure; hold on; 
% plot(mids{1}(1,:),mids{2}(1,:),'b'); 
% plot(mids{1}(2,:),mids{2}(2,:),'r');
% for i = 1:length(mids{1}) 
%     plot(mids{1}(1,i),mids{2}(1,i),'.','Color',1-[1 1 1]*i*(1/length(mids{1})),'MarkerSize',30); 
%     plot(mids{1}(2,i),mids{2}(2,i),'o','Color',1-[1 1 1]*i*(1/length(mids{1})),'MarkerSize',10); 
% end
% plot(mids{1}(1,length(loop_ranges{1})),mids{1}(2,length(loop_ranges{1})),'bo','MarkerSize',12);
% 
% plot([-5 5],[-5 5],'k--')
%%
% lev = cumsum(mids{1}(1,:))-cumsum(mids{2}(1,:));
% hev = cumsum(mids{1}(2,:))-cumsum(mids{2}(2,:));
% elev = cumsum(mids{3}(1,:))-cumsum(mids{4}(1,:));
% ehev = cumsum(mids{3}(2,:))-cumsum(mids{4}(2,:));
% 
% ilev = (mids{1}(1,:))-(mids{2}(1,:));
% ihev = (mids{1}(2,:))-(mids{2}(2,:));
% ielev = (mids{3}(1,:))-(mids{4}(1,:));
% iehev = (mids{3}(2,:))-(mids{4}(2,:));
% 
% figure; hold on; 
% plot(lev,'b.-'); 
% plot(hev,'r.-');
% 
% plot(elev,'k.-'); 
% plot(ehev,'k.-');
% 
% figure; hold on; 
% plot(lev-elev,'b.-'); 
% plot(hev-ehev,'r.-');
%% 
% cat_av = cellfun(@(x) horzcat(x{:}),av_all,'UniformOutput',0);
% cat_low = cellfun(@(x) nanmean(reshape(x(like_ind{1},:),[],1)),cat_av);
% %
% bfunc = @(x) nanmean(nanmean(x(x(:,end)==2,1))-nanmean(x(x(:,end)==1,1))); 
% [vs,is] = deal(cell(length(like_ind),1));
% [uncdiff,uncdiff_n] = deal(cell(length(voif_all),1));
% for i = 1:length(voif_all) %PD/OD/ORTH
%     for j = 1:length(voif_all{i}) %Time
%         for k = 1:length(like_ind)
%             vs{k} = reshape(voif_all{i}{j}(like_ind{k},:),[],1);
%             vs{k}(isnan(vs{k})) = [];
%             is{k} = k.*ones(length(vs{k}),1);
%         end
%         
%         totvmat = [vertcat(vs{:}), vertcat(is{:})];
%         
%         [uncdiff{i}(1,j),uncdiff{i}(2,j)] = boot_bounds(1000,bfunc,totvmat,2.5,97.5);
%         clc; fprintf('%d/%d - %d/%d\n',i,length(voif_all),j,length(voif_all{i}));
%         
%         uncdiff{i}(3,j) = bfunc(totvmat);
%         
%         uncdiff_n{i} = uncdiff{i}./cat_low(i);
%     end
% end
% %%
% plotmark = {'o-','s-','x-'};
% xsplot = [xsforplotT xsforplotG];
% figure; hold on; 
% for i = 1:length(voif_all)
% %     patch([xsplot fliplr(xsplot)],[uncdiff{i}(1,:),fliplr(uncdiff{i}(2,:))],'k',...
% %         'FaceAlpha',0.25,'EdgeAlpha',0);
% %     plot(xsplot,uncdiff{i}(3,:),plotmark{i},'Color','k');
%     
%     patch([xsplot fliplr(xsplot)],[uncdiff_n{i}(1,:),fliplr(uncdiff_n{i}(2,:))],'k',...
%         'FaceAlpha',0.25,'EdgeAlpha',0);
%     plot(xsplot,uncdiff_n{i}(3,:),plotmark{i},'Color','k');
% end
% plot(xsplot,zeros(size(xsplot)),'k--');

