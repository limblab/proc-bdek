brain_area = 'PMd';

loop_alignments = {'target','go'};
% loop_alignments = {'target'};
loop_ranges = {[-200:100:800],[-100:100:500]};
dtime = unique(diff(loop_ranges{1}));
% loop_ranges = {[50 250 600 800],[0 200]};
% loop_ranges = {[-200 0 50 250 600 800]};
boot_type = 'trials';
FR_thresh = 2;

CO_index = 1;
if exist('pdi','var')
    prediction_day_indices = pdi;
else
    prediction_day_indices = 2;
end
reachdir_col = 10;

colors_of_plot = {'k','b','r','c'};
faceA = .25;
edgeA = .25;
% Clean up
% for i = 1:length(alldays)
%     alldays(i).tt(floor(10*(alldays(i).tt(:,3)-floor(alldays(i).tt(:,3))))==9,:) = [];
%     alldays(i).tt(isnan(alldays(i).tt(:,3)),:) = [];
% end

[percdiff, xsfp, diffb, tds, totd,neurons,FD,PDneurs,ODneurs,ORTHneurs,PD,PD1,TA,...
    Firings,dfrompd,mbound,excited,suppressed] = deal(cell(1,length(loop_alignments)));
% Initialize
co_units = eval(sprintf('alldays(1).%s_units',brain_area));
if strcmp(brain_area,'M1'); co_units(end-1:end) =[]; end
baseraster = cell(length(co_units),1);
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
        
                [cents,tune,tune_low,tune_high] = co_tuning(co_units,[nan; alldays(CO_index).tt(1:end-1,10)],...
                alldays(CO_index).tt,t1,t2,tune_align,i);
%                 [wrapped_cents,wrapped_tune,p,wrapped_low,wrapped_high] = co_tuning_VM(co_units,[nan; alldays(CO_index).tt(1:end-1,reachdir_col)],...
%                 alldays(CO_index).tt,t1,t2,tune_align,i,0);
               
                baseraster{i,:}{ranges} = raster_plot(co_units{i},alldays(CO_index).tt,[t1 t2]./1000,tune_align,0,'none');
                
            else
                
                [cents,tune,tune_low,tune_high] = co_tuning(co_units,alldays(CO_index).tt(:,10),...
                alldays(CO_index).tt,t1,t2,tune_align,i);
            
%                 [wrapped_cents,wrapped_tune,p,wrapped_low,wrapped_high] = co_tuning_VM(co_units,alldays(CO_index).tt(:,reachdir_col),...
%                 alldays(CO_index).tt,t1,t2,tune_align,i,0);

    %             [cents,tune,tune_low,tune_high] = co_tuning(co_units,alldays(CO_index).tt(:,2),...
    %             alldays(CO_index).tt,0,200,tune_align,i);
            end
            if length(wrapped_low)~=length(wrapped_tune); wrapped_low=wrapped_tune; wrapped_high=wrapped_tune; end
            % Smooth and fix edges on tuning curves
%             front_cent = cents;
%             back_cent = cents;
% 
%             front_tune = tune; front_low = tune_low; front_high = tune_high;
%             back_tune = tune; back_low = tune_low; back_high = tune_high;
% 
%             tune_cents = [back_cent-2*pi cents front_cent+2*pi];
%             tune_pad = [back_tune;tune;front_tune];
%             low_pad = [back_low;tune_low;front_low]; 
%             high_pad = [back_high;tune_high;front_high];
% 
%             interp_cents = tune_cents(1):.01:tune_cents(end);
%             tune_interp = interp1(tune_cents,tune_pad,interp_cents);
%             low_interp = interp1(tune_cents,low_pad,interp_cents);
%             high_interp = interp1(tune_cents,high_pad,interp_cents);
% 
%             smooth_tune = smooth(tune_interp,100);
%             smooth_low = smooth(low_interp,100);
%             smooth_high= smooth(high_interp,100);
% 
%             wrapped_cents = interp_cents(interp_cents >=0 & interp_cents< 2*pi);
%             wrapped_tune = smooth_tune(interp_cents >= 0 & interp_cents <2*pi); 
%             wrapped_low = smooth_low(interp_cents >= 0 & interp_cents <2*pi); 
%             wrapped_high = smooth_high(interp_cents >= 0 & interp_cents <2*pi);
% 
%             neurons{loopthrough}{ranges}{i}.tuning = wrapped_tune;
%             neurons{loopthrough}{ranges}{i}.tuning_low = wrapped_low;
%             neurons{loopthrough}{ranges}{i}.tuning_high = wrapped_high;
%             
%             PD{loopthrough}{ranges}{i} = wrapped_cents(wrapped_tune==max(wrapped_tune));
%             if length(PD{loopthrough}{ranges}{i}) < 2; PD1{loopthrough}{ranges}(i,:) = PD{loopthrough}{ranges}{i}; else
%             PD1{loopthrough}{ranges}(i,:) = NaN; end
        
            PD = wrapped_cents(wrapped_tune==max(wrapped_tune));
            
            if length(PD)>1 || isempty(PD)
                PD1{loopthrough}{ranges}(i,:) = NaN;%
            else
                PD1{loopthrough}{ranges}(i,:) = PD;
            end

            neurons{loopthrough}{ranges}{i}.tuning = wrapped_tune';
            neurons{loopthrough}{ranges}{i}.tuning_low = wrapped_low';
            neurons{loopthrough}{ranges}{i}.tuning_high = wrapped_high';
        end
        
    reppd = repmat(best_PDS',size(alldays(prediction_day_indices).tt,1),1);
    %reppd = repmat(PD1{loopthrough}{ranges}',size(alldays(prediction_day_indices).tt,1),1);
    reprch = repmat(alldays(prediction_day_indices).tt(:,10),1,size(reppd,2));
    
    dfrompd{loopthrough}{ranges} = circ_dist(reppd,reprch);
    
    PDinds = abs(circ_dist(reppd,reprch)) <= pi/4; 
    ODinds  = abs(circ_dist(reppd,reprch)) >= 3*pi/4; 
    ORTHinds = abs(circ_dist(reppd,reprch))>pi/4 & abs(circ_dist(reppd,reprch))<3*pi/4; 
     
    PDneurs{loopthrough}{ranges} = nan(size(PDinds)); PDneurs{loopthrough}{ranges}(PDinds)=1;
    ODneurs{loopthrough}{ranges} = nan(size(ODinds)); ODneurs{loopthrough}{ranges}(ODinds)=1;
    ORTHneurs{loopthrough}{ranges} = nan(size(ORTHinds)); ORTHneurs{loopthrough}{ranges}(ORTHinds)=1;
        
    end
    
end
static_PD = PDneurs{1}{end};
static_OD = ODneurs{1}{end};
static_ORTH = ORTHneurs{1}{end};

allneurs = horzcat(neurons{:});
[time_tuning_mat,tunearrays] = deal(cell(1,length(allneurs)));
for i = 1:length(allneurs)
    time_tunings = cellfun(@(x) x.tuning',allneurs{i},'UniformOutput',0);
    tunearrays{i} = cell2mat(cellfun(@(x) x.tuning,allneurs{i},'UniformOutput',0))';
    time_tuning_mat{i} = vertcat(time_tunings{:});
end
alltimetunes = horzcat(time_tuning_mat{:});

baselinematrix = cellfun(@(x) cell2mat(horzcat(x{:})),baseraster,'UniformOutput',0);
baseline_levels = cell2mat(cellfun(@(x) nanmean(nansum(x,2)./(size(x,2)/1000),1),baselinematrix,'UniformOutput',0));

%%
figure; hold on;
for loopthrough = 1:length(loop_alignments)
    
    tune_align = loop_alignments{loopthrough};
    tune_ranges = loop_ranges{loopthrough};
    tds{loopthrough} = diff(tune_ranges);
    time_bins = tune_ranges;
    time_align = tune_align;

    unit_counts = cell(length(prediction_day_indices),1);
    for z = 1:length(prediction_day_indices)

        prediction_day = prediction_day_indices(z);
        day_units = eval(sprintf('alldays(1).%s_units',brain_area));
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

            unit_counts{z}{q} = rast./repmat(diff(time_bins)/1000,size(rast,1),1); % in FIRING RATE
        end
    end

    %
    [expect_counts,firing_perc, trial_counts, lowerbound, upperbound, ...
        midbound, all_diffs, firing_diffs, good_neurso, good_neurs, ...
        firing_absolute,straight_difference] = deal(cell(length(prediction_day_indices),1));
    
    for bin = 1:length(time_bins) - 1

        tune_range = tune_ranges(bin:bin+1);
        t1 = tune_range(1);
        t2 = tune_range(2);

        clc; fprintf('bin %d/%d\n',bin,length(time_bins)-1);
        unit_cell = struct2cell(vertcat(neurons{loopthrough}{bin}{:}));    
        tuning_array = horzcat(unit_cell{1,:})'; % in FIRING RATE
        
        TA{loopthrough}{bin} = tuning_array;

        modulations = (max(tuning_array,[],2) - min(tuning_array,[],2))'; % in FR
%         good_neurso{bin} = find(modulations > FR_thresh);
        good_neurso{bin} = find(modulations > FR_thresh);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %good_neurs{bin} = good_neurso{bin}(ismember(good_neurso{bin},far_in));
        good_neurs{bin} = good_neurso{bin};
        %good_neurs{bin} = 119;
        %good_neurs{bin} = good_neurso{bin}(ismember(good_neurso{bin},GOOD_list)); 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        gns = nan(1,length(modulations)); gns(good_neurs{bin}) = 1;

        
        %colors_of_plot = {'c','g','m','r'};
        moduls = cell(length(prediction_day_indices),1);
        FR_change = moduls;
        like_ind = moduls;
        for z = 1:length(prediction_day_indices)

            day_pred = prediction_day_indices(z);
            tt_pred = alldays(day_pred).tt;  
            sl_pred = alldays(day_pred).slices;

            moduls{z} = repmat(modulations,size(tt_pred,1),1);
      
            day_units = eval(sprintf('alldays(1).%s_units',brain_area));
            if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end

            for i = 1:size(tt_pred,1)   
                for q = 1:length(day_units)       

                    trial_counts{z}{bin}(i,q) = unit_counts{z}{q}(i,bin); % in FIRING RATE
                    %move_aligned_thetas = abs(circ_dist(wrapped_cents,tt_pred(i,9)));
                    move_aligned_thetas = abs(circ_dist(wrapped_cents,tt_pred(i,10)));
                    move_theta_index = find(move_aligned_thetas==min(move_aligned_thetas),1,'first');
                    expect_counts{z}{bin}(i,q) = tuning_array(q,move_theta_index); % in FIRING RATE
                    
                    %%% RECRUIT BY SLICE
% %                     slice_count = zeros(1,size(sl_pred,2));
% %                     for sls = 1:size(sl_pred,2)
% %                         move_aligned_thetas = abs(circ_dist(wrapped_cents,sl_pred(i,sls)));
% %                         move_theta_index = find(move_aligned_thetas==min(move_aligned_thetas),1,'first');
% %                         slice_count(sls) = tuning_array(q,move_theta_index);
% %                     end
% %                    % expect_counts{z}{bin}(i,q) = sum(slice_count)/5; % in FIRING RATE 
% %                     expect_counts{z}{bin}(i,q) = max(slice_count); % in FIRING RATE
                    %%
                end
            end
            
            if loopthrough ==1 && loop_ranges{1}(bin)<0
                excited{loopthrough}{bin} = ones(size(expect_counts{z}{bin}));
                suppressed{loopthrough}{bin} = ones(size(expect_counts{z}{bin}));
            else
                excited{loopthrough}{bin} = double(expect_counts{z}{bin} - repmat(baseline_levels',length(tt_pred),1) > 0); 
                excited{loopthrough}{bin}(excited{loopthrough}{bin}==0)=NaN;
                suppressed{loopthrough}{bin} = double(expect_counts{z}{bin} - repmat(baseline_levels',length(tt_pred),1) < 0); 
                suppressed{loopthrough}{bin}(suppressed{loopthrough}{bin}==0)=NaN;
            end
            
            
            straight_difference{z}{bin} = (trial_counts{z}{bin}(:,good_neurs{bin}) - ...
                expect_counts{z}{bin}(:,good_neurs{bin}));
            all_diffs{z}{bin} = (trial_counts{z}{bin}-expect_counts{z}{bin})./moduls{z};
            firing_diffs{z}{bin} = straight_difference{z}{bin}./moduls{z}(:,good_neurs{bin});
            firing_perc{z}{bin} = (trial_counts{z}{bin}(:,good_neurs{bin}) - ...
                repmat(min(tuning_array(good_neurs{bin},:),[],2)',length(tt_pred),1))./moduls{z}(:,good_neurs{bin});
            firing_absolute{z}{bin} = trial_counts{z}{bin}(:,good_neurs{bin});
            straight_difference{z}{bin} = (trial_counts{z}{bin}(:,good_neurs{bin}) - ...
                expect_counts{z}{bin}(:,good_neurs{bin}));

            %operating = (trial_counts{z}{bin}-repmat(operating_mins',length(tt_pred),1))./repmat(operating_ranges',length(tt_pred),1);
            %operating_valid = operating(:,good_neurs{bin});
            
            proportion_excited = double((trial_counts{z}{bin} - repmat(baseline_levels',length(tt_pred),1))>0);
%             
            pe = proportion_excited; pe(pe==0)=NaN;
%             ECS{loopthrough}{bin} = expect_counts{z}{bin};
%             TCS{loopthrough}{bin} = trial_counts{z}{bin};

            from_baseline = trial_counts{z}{bin}-repmat(baseline_levels',length(tt_pred),1);
%             from_baseline = 
            pe_frombaseline = from_baseline.*pe;
            
%             FR_change{z} = nanmean((trial_counts{z}{bin}(:,good_neurs{bin}) - ...
%                 expect_counts{z}{bin}(:,good_neurs{bin}))./repmat(operating_ranges(good_neurs{bin})',length(tt_pred),1),2);
%             FR_change{z} = nanmean(straight_difference{z}{bin}.*suppressed{loopthrough}{bin}(:,good_neurs{bin}),2);
%             FR_change{z} = nanmean(proportion_excited,2);
            FR_change{z} = nanmean(from_baseline.*static_PD,2);%.*static_PD,2); 
            %FR_change{z} = var(trial_counts{z}{bin},[],2)-var(expect_counts{z}{bin},[],2);
            %FR_change{z} = nanmean(firing_diffs{z}{bin},2);
            %FR_change{z} = nanmean(trial_counts{z}{bin}(:,good_neurs{bin}) - expect_counts{z}{bin}(:,good_neurs{bin}),2);
            %FR_change{z} = nanmean(operating_valid,2);
            %FR_change{z} = nanvar(operating_valid,[],2);
          
            %FR_change{z} = nansum(firing_absolute{z}{bin}.*diff(time_bins(1:2)/1000),2);

            %FR_change{z} = var(firing_diffs{z}{bin},[],2);%./mean(firing_diffs{z}{bin},2);
            %FR_change{z} = nanmean(firing_absolute{z}{bin},2);
            %FR_change{z} = mean(firing_perc{z}{bin},2);
            %FR_change{z} = mean(firing_absolute{z}{bin},2);
            
            gnsrep = repmat(gns,size(alldays(prediction_day_indices(z)).tt,1),1);
            gnsrep(gnsrep==0)=nan;
%             FDnear = gnsrep.*static_near.*all_diffs{z}{bin};
%             FDfar = gnsrep.*static_far.*all_diffs{z}{bin};

            %FR_change{z} = firing_absolute{z}{bin};
            Firings{loopthrough}{bin} = trial_counts{z}{bin};
            FRDS{loopthrough}{bin} = all_diffs{z}{bin};
            %FR_change{z} = nanmean(FDfar,2);

            set_of_inds = unique(alldays(prediction_day_indices(z)).tt(:,3));
            set_of_inds = flipud(set_of_inds);

            for qq = 1:length(set_of_inds)

                if t1 <= 0 && t2 <= 0 && loopthrough == 1
                    %like_ind{z}{qq} = find(tt_pred(:,3)==set_of_inds(qq));
                    like_ind{z}{qq} = find(ismember(tt_pred(:,3),set_of_inds));
                else
                    like_ind{z}{qq} = find(ismember(tt_pred(:,3),set_of_inds(qq)));
                end

                midbound{z}(qq,bin) = mean(FR_change{z}(like_ind{z}{qq}));  
                FD{loopthrough}{z}{qq}(:,bin) = nanmean(all_diffs{z}{bin}(like_ind{z}{qq},:))';
                mbound{loopthrough}(qq,bin) = midbound{z}(qq,bin);

                if strcmp(boot_type,'trials')

                    [lowerbound{z}(qq,bin),upperbound{z}(qq,bin)] = ...
                        boot_bounds(10000,@mean,FR_change{z}(like_ind{z}{qq}),2.5,97.5);     

                elseif strcmp(boot_type,'neurons')

                    boot_across_neurs = mean(firing_diffs{z}{bin}(like_ind{z}{qq},:),1);

                    [lowerbound{z}(qq,bin),upperbound{z}(qq,bin)] = ...
                        boot_bounds(10000,@mean,boot_across_neurs,2.5,97.5);    

                elseif strcmp(boot_type,'all')

                    neuron_trials = firing_diffs{z}{bin}(like_ind{z}{qq},:);
                    cat_neur_tri = reshape(neuron_trials,1,numel(neuron_trials));

                    [lowerbound{z}(qq,bin),upperbound{z}(qq,bin)] = ...
                        boot_bounds(10000,@mean,cat_neur_tri,2.5,97.5);     
                else

                    fprintf('Improper bootstrap dimension: Bootstrapping across all\n'); 

                    neuron_trials = firing_diffs{z}{bin}(like_ind{z}{qq},:);
                    cat_neur_tri = reshape(neuron_trials,1,numel(neuron_trials));

                    [lowerbound{z}(qq,bin),upperbound{z}(qq,bin)] = ...
                        boot_bounds(10000,@mean,cat_neur_tri,2.5,97.5);   

                end

                totd{loopthrough}{qq}(:,bin) = [midbound{z}(qq,bin); lowerbound{z}(qq,bin); upperbound{z}(qq,bin)];
              
            end
        end
    end
    %
    dayind = 1;

    if loopthrough == 1

        xsforplotT = .5*(time_bins(1:end-1)+time_bins(2:end));
        max_targ_x = max(xsforplotT);
        for i = 1:length(unique(alldays(prediction_day_indices).tt(:,3)))
            plot(xsforplotT,midbound{dayind}(i,:),colors_of_plot{i});
            patch([xsforplotT fliplr(xsforplotT)],[lowerbound{dayind}(i,:) fliplr(upperbound{dayind}(i,:))],...
                colors_of_plot{i},'FaceAlpha',faceA,'EdgeAlpha',edgeA);
        end
        xsfp{1} = xsforplotT;

    else

        xsatgo = .5*(time_bins(1:end-1)+time_bins(2:end));
        xsforplotG = xsatgo + max_targ_x + 100;
        min_go_x = min(xsatgo);

        for i = 1:length(unique(alldays(prediction_day_indices).tt(:,3)))
            plot(xsforplotG,midbound{dayind}(i,:),colors_of_plot{i});
            patch([xsforplotG fliplr(xsforplotG)],[lowerbound{dayind}(i,:) fliplr(upperbound{dayind}(i,:))],...
                colors_of_plot{i},'FaceAlpha',faceA,'EdgeAlpha',edgeA);
        end 
        xsfp{2} = xsforplotG;
    end
end
%%
FDc = horzcat(FD{:});
[pos_totd] = deal(cell(1,2));
for prepost = 1:length(loop_alignments)
    positiveinds = xsfp{prepost} > 0;
    pos_totd{prepost} = cellfun(@(x) x(:,positiveinds), totd{prepost},'UniformOutput',0);
end
% [integs,allint,totint,nbn_integs] = deal(cell(length(unique(alldays(dayind+1).tt(:,3))),1));
% for lik = 1:length(unique(alldays(dayind+1).tt(:,3)))
%     integs{lik} = cellfun(@(x) trapz(x{lik},2),pos_totd,'UniformOutput',0);
%     nbn_integs{lik} = cellfun(@(x) trapz(x{lik},2),FDc,'UniformOutput',0);
%     
%     allint{lik} = sum(horzcat(integs{lik}{:}),2);
%     totint{lik} = sum(horzcat(nbn_integs{lik}{:}),2);
% end
% errorbar_integral = horzcat(allint{:}).*(unique(diff(xsforplotT))./1000);
% nbn_integral = horzcat(totint{:}).*(unique(diff(xsforplotT))./1000);
% 
% catted = vertcat(FDc{:});
% [TI,total_int] = deal(cell(1,size(catted,2)));
% for i = 1:size(catted,2)
%     TI{i} = horzcat(catted{:,i});
%     integs{i} = TI{i};
%     integs{i}(isinf(integs{i})) = NaN;
%     total_int{i} = nansum(integs{i},2);
% end

%ylim([-.5 1.5]);
plot(max_targ_x*[1 1],[-.5 1.5],'k--');
plot([0 0],[-.5 1.5],'k--');
plot((max_targ_x+min_go_x+100)*[1 1],[-.5 1.5],'k--');
xlabel('Time From Target On (ms)','FontSize',14);
ylabel('Change in FR (% modulation)','FontSize',14);
title(sprintf('%s (%s) - %s/%s/%s',monkey,brain_area,MO,DA,YE),'FontSize',16);

% figure; hold on;
% bar(xsforplotT,percdiff{1}); bar(xsforplotG,percdiff{2});
% plot(max_targ_x*[1 1],[0 1],'k--');
% plot([0 0],[0 1],'k--');
% xlabel('bin'); ylabel('% Different (H - L)'); 

if 0
%%

neurons_all = [neurons{1} neurons{2}];

[tunings, deviations] = deal(cell(1,length(neurons_all)));
moddepth = zeros(length(neurons_all{1}),length(tunings));
for i = 1:length(neurons_all); 
    for j = 1:length(neurons_all{1}); 
        tunings{i}(j,:) = neurons_all{i}{j}.tuning'; 
    end
end
basemat = horzcat(tunings{1:4});
baselevels = nanmean(basemat,2);
baselevel_mat = repmat(baselevels,1,628);
for i = 1:length(tunings); 
    deviations{i} = tunings{i}-baselevel_mat; 
    moddepth(:,i) = max(tunings{i},[],2)-min(tunings{i},[],2);
end
%
[mdl,mdh] = boot_bounds(1000,@nanmean,moddepth,2.5,97.5);

x{1} = 0.5*(loop_ranges{1}(1:end-1)+loop_ranges{1}(2:end));
x{2} = 0.5*(loop_ranges{2}(1:end-1)+loop_ranges{2}(2:end))+1000;
figure; hold on; plot([x{1} x{2}],nanmean(moddepth),'k'); 
patch([[x{1} x{2}] fliplr([x{1} x{2}])],[mdl' fliplr(mdh')],'k','FaceAlpha',0.25,'EdgeAlpha',0);
end

%%
% figure; hold on; 
% plot([low_excite(2,1) up_excite(2,1)],[1 1]*mid_suppress(2,1),'b','LineWidth',2); plot([1 1]*mid_excite(2,1),[low_suppress(2,1),up_suppress(2,1)],'b','LineWidth',2)
% plot([low_excite(2,2) up_excite(2,2)],[1 1]*mid_suppress(2,2),'b'); plot([1 1]*mid_excite(2,2),[low_suppress(2,2),up_suppress(2,2)],'b')
% plot([low_excite(3,1) up_excite(3,1)],[1 1]*mid_suppress(3,1),'r','LineWidth',2); plot([1 1]*mid_excite(3,1),[low_suppress(3,1),up_suppress(3,1)],'r','LineWidth',2)
% plot([low_excite(3,2) up_excite(3,2)],[1 1]*mid_suppress(3,2),'r'); plot([1 1]*mid_excite(3,2),[low_suppress(3,2),up_suppress(3,2)],'r')
% plot(mid_excite(2,:),mid_suppress(2,:),'k'); plot(mid_excite(3,:),mid_suppress(3,:),'k')

