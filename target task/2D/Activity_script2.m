% for i = 1:length(alldays)
%     alldays(i).tt(floor(10*(alldays(i).tt(:,3)-floor(alldays(i).tt(:,3))))==9,:) = [];
%     alldays(i).tt(isnan(alldays(i).tt(:,3)),:) = [];
% end
tune_ranges = [50 250];
tune_align = 'go';
brain_area = 'PMd_narrow';
MODthresh = 0;
if exist('pdi','var')
    prediction_day_indices = pdi;
else
    prediction_day_indices = 2;
end
CO_index = 1;
group_col = 3; 
reachdir_col = 10;

[neurons, PD,PD1,nears,fars,nearneurs,farneurs,DFPD,VM_P] = deal(cell(length(tune_ranges)-1,1));
for ranges = 1:length(tune_ranges)-1
    tune_range = tune_ranges(ranges:ranges+1);
    % Initialize
    co_units = eval(sprintf('alldays(1).%s_units',brain_area));
    if strcmp(brain_area,'M1'); co_units(end-1:end) =[]; end
    %neurons = cell(length(co_units),length(tune_range)-1);

    for i = 1:length(co_units) % Loop through units
         clc; fprintf('%d/%d Tuning: (%d/%d)\n',ranges,length(tune_ranges)-1,i,length(co_units));

        % set time bin edges
        t1 = tune_range(1);
        t2 = tune_range(2);

        [cents,tune,tune_low,tune_high] = co_tuning(co_units,alldays(CO_index).tt(:,reachdir_col),...
            alldays(CO_index).tt,t1,t2,tune_align,i);
       
        %%Smooth and fix edges on tuning curves
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

        neurons{ranges}{i}.tuning = wrapped_tune;
        neurons{ranges}{i}.tuning_low = wrapped_low;
        neurons{ranges}{i}.tuning_high = wrapped_high;
      
        PD{ranges}{i} = wrapped_cents(wrapped_tune==max(wrapped_tune));
        if length(PD{ranges}{i}) < 2; PD1{ranges}(i,:) = PD{ranges}{i}; else
            PD1{ranges}(i,:) = NaN; end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         [wrapped_cents,wrapped_tune,p,wrapped_low,wrapped_high] = co_tuning_VM(co_units,alldays(CO_index).tt(:,reachdir_col),...
%             alldays(CO_index).tt,t1,t2,tune_align,i,1);
%         
%         VM_P{ranges}{i} = p;
%         
%         PD{ranges}{i} = wrapped_cents(wrapped_tune==max(wrapped_tune));
%         PD{ranges}{i} = wrapped_cents(abs(circ_dist(wrapped_cents,p(end)))==min(abs(circ_dist(wrapped_cents,p(end)))));
%         PD1{ranges}(i,:) = PD{ranges}{i};
%         if length(PD{ranges}{i})>1
%             PD1{ranges}(i,:) = NaN;%PD{ranges}{1};
%         else
%             PD1{ranges}(i,:) = PD{ranges}{i};
%         end
%         
%         neurons{ranges}{i}.tuning = wrapped_tune;
%         neurons{ranges}{i}.tuning_low = wrapped_low;
%         neurons{ranges}{i}.tuning_high = wrapped_high;
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    end
%     nears{ranges} = find(abs(circ_dist(PD1{ranges},circ_mean(alldays(prediction_day_indices).tt(:,reachdir_col))))< 2*pi/4);
%     fars{ranges} = find(abs(circ_dist(PD1{ranges},circ_mean(alldays(prediction_day_indices).tt(:,reachdir_col))))> 2*pi/4);
%     
    reppd = repmat(PD1{ranges}',size(alldays(prediction_day_indices).tt,1),1);
    %reprch = repmat(alldays(prediction_day_indices).tt(:,reachdir_col),1,size(reppd,2));
    reprch = pi/2*ones(size(reppd));
    
    dfrompd = circ_dist(reppd,reprch);
    
%     nearinds = abs(dfrompd) < pi/2; 
%     farinds  = abs(dfrompd) > pi/2; 
%      
%     nearneurs{ranges} = nan(size(nearinds)); nearneurs{ranges}(nearinds)=1;
%     farneurs{ranges} = nan(size(farinds)); farneurs{ranges}(farinds)=1;
    
    DFPD{ranges} = dfrompd;
    
end
centers = round(wrapped_cents*1000)./1000;


%%
figure; hold on; 
time_bins = tune_ranges;
time_align = tune_align;

movemean = circ_mean(alldays(prediction_day_indices).tt(:,reachdir_col));
xsforplot = .5*(time_bins(1:end-1)+time_bins(2:end));
 
[unit_counts,mintrial,maxtrial] = deal(cell(length(prediction_day_indices),1));
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

        unit_counts{z}{q} = rast./repmat(diff(time_bins)/1000,size(rast,1),1);
        mintrial{z}{q} = min(unit_counts{z}{q},[],2);
        maxtrial{z}{q} = max(unit_counts{z}{q},[],2);
    end
end
allmin = horzcat(mintrial{z}{:});
allmax = horzcat(maxtrial{z}{:});

percnegfunc = @(F,mint,maxt) 2*F./(maxt-mint) - (maxt+mint)./(maxt-mint);
%%

[expect_counts, firing_perc, trial_counts, lowerbound, upperbound, ...
    midbound, firing_diffs, good_neurs, dfr,closemove,farmove, efiring_perc,...
    firing_absolute, firing_minmax, firing_expect, midbounde, TA, MODS,...
    spike_count_raw,neur_fd,modulation_bin,FD_allneurs,firing_percneg,...
    efiring_percneg,goods,E_perc,I_perc,ex_neurs,in_neurs] = deal(cell(length(prediction_day_indices),1));

excites = cell(length(time_bins)-1,1);
for bin = 1:length(time_bins) - 1
    clc; fprintf('bin %d/%d\n',bin,length(time_bins)-1);
    unit_cell = struct2cell(vertcat(neurons{bin}{:}));    
    tuning_array = horzcat(unit_cell{1,:})';
    %tuning_array = vertcat(unit_cell{1,:});

    modulations = (max(tuning_array,[],2) - min(tuning_array,[],2))';
    good_neurs{bin} = find(modulations > MODthresh);
    gns = nan(1,length(modulations)); gns(good_neurs{bin}) = 1;
    goods{bin} = repmat(gns,size(eval(sprintf('alldays(%d).tt',prediction_day_indices)),1),1);
    TA{bin} = tuning_array;
    baselines = nanmean(TA{1},2);
    
    if (time_bins(bin)+time_bins(bin+1)) <= 0
        excites{bin} = 2*ones(size(TA{bin},1),1);
    else
        excites{bin} = nanmean(TA{bin},2)>baselines;
    end
    
    modulation_bin{1}{bin} = modulations; 
    
    [closemove{bin},farmove{bin}] = deal(zeros(size(tt_pred,1),length(day_units)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %good_neurs{bin} = good_neurs{bin}(ismember(good_neurs{bin},fars{end}));
    %good_neurs{bin} = good_neurs{bin}(ismember(good_neurs{bin},GOOD_list));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %colors_of_plot = {'c','r','r','c'};
    colors_of_plot = {'b','c','r','m'};
    %colors_of_plot = {'c','b','m','r'};
    %figure; hold on;
    [moduls, FR_change, like_ind, FR_expect] = deal(cell(length(prediction_day_indices),1));
    for z = 1:length(prediction_day_indices)
        
        day_pred = prediction_day_indices(z);
        tt_pred = alldays(day_pred).tt;  
       
        moduls{z} = repmat(modulations,size(tt_pred,1),1);
        MODS{bin} = moduls{z};

        day_units = eval(sprintf('alldays(1).%s_units',brain_area));
        if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end
 
        
        move_theta_index = zeros(size(tt_pred,1),1);
        for i = 1:size(tt_pred,1)   
            
            closemove{bin}(i,:) = abs(circ_dist(PD1{bin},alldays(prediction_day_indices).tt(i,reachdir_col))) < pi/2;
            farmove{bin}(i,:) = abs(circ_dist(PD1{bin},alldays(prediction_day_indices).tt(i,reachdir_col))) > pi/2;
            closemove{bin}(:,~ismember(1:length(day_units),good_neurs{bin})) = 0;
            farmove{bin}(:,~ismember(1:length(day_units),good_neurs{bin})) = 0;
            
            %move_aligned_thetas = abs(circ_dist(wrapped_cents,tt_pred(i,reachdir_col)+pi));
            move_aligned_thetas = abs(circ_dist(wrapped_cents,tt_pred(i,reachdir_col)));
            move_aligned_thetas2 = abs(circ_dist(wrapped_cents,tt_pred(i,9)));
            move_theta_index(i) = find(move_aligned_thetas==min(move_aligned_thetas),1,'first');
            %move_theta_index2(i) = find(move_aligned_thetas2==min(move_aligned_thetas2),1,'first');
            
            for q = 1:length(day_units)       
                trial_counts{z}{bin}(i,q) = unit_counts{z}{q}(i,bin);
                expect_counts{z}{bin}(i,q) = tuning_array(q,move_theta_index(i));
                %expect_counts2{z}{bin}(i,q) = tuning_array(q,move_theta_index2(i));
            end
            
        end
%       
        expect_counts{z}{bin}(expect_counts{z}{bin}==0)= nan;
        %expect_counts2{z}{bin}(expect_counts2{z}{bin}==1)=nan;
        
%         likelihood_r{z}{bin} = log(poisspdf(trial_counts{z}{bin},expect_counts{z}{bin})./...
%                                poisspdf(trial_counts{z}{bin},expect_counts2{z}{bin}));
        
        firing_diffs{z}{bin} = (trial_counts{z}{bin}(:,good_neurs{bin}) - ...
            expect_counts{z}{bin}(:,good_neurs{bin}))./moduls{z}(:,good_neurs{bin});
%         firing_diffsNN{z}{bin} = trial_counts{z}{bin}.*goods{bin} - ...
%             expect_counts{z}{bin}.*goods{bin};
%             
%         firing_diffs{z}{bin} = (trial_counts{z}{bin}.*goods{bin} - ...
%             expect_counts{z}{bin}.*goods{bin})./(moduls{z}.*goods{bin});       
%         
%         firing_perc{z}{bin} = (trial_counts{z}{bin}(:,good_neurs{bin}) - ...
%             repmat(min(tuning_array(good_neurs{bin},:),[],2)',length(tt_pred),1))./moduls{z}(:,good_neurs{bin});
%         
%         firing_percall{z}{bin} = (trial_counts{z}{bin}.*goods{bin} - ...
%             repmat(min(tuning_array,[],2)',length(tt_pred),1))./(moduls{z}.*goods{bin});
%         
%         efiring_perc{z}{bin} = (expect_counts{z}{bin}(:,good_neurs{bin}) - ...
%             repmat(min(tuning_array(good_neurs{bin},:),[],2)',length(tt_pred),1))./moduls{z}(:,good_neurs{bin});
%         firing_absolute{z}{bin} = trial_counts{z}{bin}(:,good_neurs{bin});
%         spike_count_raw{z}{bin} = trial_counts{z}{bin}(:,good_neurs{bin})*(diff(time_bins(1:2))/1000);
%         firing_minmax{z}{bin} = (trial_counts{z}{bin}(:,good_neurs{bin}) - ...
%             allmin(:,good_neurs{bin}))./(allmax(:,good_neurs{bin}) - allmin(:,good_neurs{bin}));
%         firing_expect{z}{bin} = expect_counts{z}{bin}(:,good_neurs{bin});
%         
%         firing_percneg{z}{bin} = percnegfunc(trial_counts{z}{bin}(:,good_neurs{bin})...
%                                             ,repmat(min(tuning_array(good_neurs{bin},:),[],2)',length(tt_pred),1)...
%                                             ,repmat(max(tuning_array(good_neurs{bin},:),[],2)',length(tt_pred),1));
%         efiring_percneg{z}{bin} = percnegfunc(expect_counts{z}{bin}(:,good_neurs{bin})...
%                                             ,repmat(min(tuning_array(good_neurs{bin},:),[],2)',length(tt_pred),1)...
%                                             ,repmat(max(tuning_array(good_neurs{bin},:),[],2)',length(tt_pred),1));
%                                         
%                     
%                                         
%         ex_neurs{bin} = good_neurs{bin}(ismember(excites{bin}(good_neurs{bin}),[1 2]));
%         in_neurs{bin} = good_neurs{bin}(ismember(excites{bin}(good_neurs{bin}),[0 2]));
%         
%         E_perc{z}{bin} = (trial_counts{z}{bin}(:,ex_neurs{bin}) - ...
%             repmat(min(tuning_array(ex_neurs{bin},:),[],2)',length(tt_pred),1))./moduls{z}(:,ex_neurs{bin});                           
%         I_perc{z}{bin} = (trial_counts{z}{bin}(:,in_neurs{bin}) - ...
%             repmat(min(tuning_array(in_neurs{bin},:),[],2)',length(tt_pred),1))./moduls{z}(:,in_neurs{bin}); 
%         
%         FD_allneurs{z}{bin} = (trial_counts{z}{bin} - expect_counts{z}{bin})./moduls{z};
        
        FR_change{z} = nanmean(firing_diffs{z}{bin},2);
        %FR_change{z} = mean(firing_perc{z}{bin},2);
        %FR_change{z} = nanmean(firing_minmax{z}{bin},2);
        %FR_change{z} = nansum(firing_absolute{z}{bin},2);
        %FR_expect{z} = nansum(firing_expect{z}{bin},2);
        %dfr{z}(:,bin) = FR_change{z};
        
%         gnsrep = repmat(gns,size(alldays(prediction_day_indices(z)).tt,1),1);
%         gnsrep(gnsrep==0)=nan;
%         FDnear = gnsrep.*nearneurs{1}.*firing_diffs{1}{bin};
%         FDfar = gnsrep.*farneurs{1}.*firing_diffs{1}{bin};
        
        %FR_change{z} = nanmean(FDnear,2);
       
        set_of_inds = unique(alldays(prediction_day_indices(z)).tt(:,group_col));
        set_of_inds = flipud(set_of_inds);
        
        %like_ind = li2;
        for qq = 1:length(set_of_inds)
            like_ind{z}{qq} = find(tt_pred(:,group_col)==set_of_inds(qq));
            %[lk{z}{qq},uk{z}{qq}] = boot_bounds(10000,@mean,FR_change{z}(like_ind{z}{qq}),2.5,97.5);   
            [lowerbound{z}(qq,bin),upperbound{z}(qq,bin)] = ...
                boot_bounds(10000,@mean,FR_change{z}(like_ind{z}{qq}),2.5,97.5);     
            midbound{z}(qq,bin) = mean(FR_change{z}(like_ind{z}{qq}));         
            
            neur_fd{z}{qq} = nanmean(firing_diffs{z}{bin}(like_ind{z}{qq},:),1);
            %midbounde{z}(qq,bin) = mean(FR_expect{z}(like_ind{z}{qq}));
        end
    end
end
%fdif = FR_diff_NbN(tt_pred,firing_diffs{1},day_units,good_neurs);
% good_neurso = good_neurs;
% significant_diffs_in_FR;
% figure; hold on; bar(xsforplot,perc_d_inbin);
% xlabel('bin'); ylabel('% Different (H - L)'); 
%
%figure;
hold on;
%subplot(1,2,2); hold on;
for i = 1:length(unique(alldays(prediction_day_indices).tt(:,group_col)))
    
    plot(xsforplot,midbound{1}(i,:),colors_of_plot{i});
    patch([xsforplot fliplr(xsforplot)],[lowerbound{1}(i,:) fliplr(upperbound{1}(i,:))],...
        colors_of_plot{i},'FaceAlpha',0.75,'EdgeAlpha',0.75);
%     
%     plot(xsforplot,midbound{dayind}(i,:) - midbounde{dayind}(i,:),colors_of_plot{i});
%     patch([xsforplot fliplr(xsforplot)],[lowerbound{dayind}(i,:) fliplr(upperbound{dayind}(i,:))] - [midbounde{dayind}(i,:) fliplr(midbounde{dayind}(i,:))],...
%         colors_of_plot{i},'FaceAlpha',0.5,'EdgeAlpha',0.5);
end 
xlabel('Time From Target On (ms)','FontSize',14);
ylabel('Change in FR (% modulation)','FontSize',14);

%ylim([-0.5 1.5]);
%%
% 
% DIFFP = cell(3,1);
% figure; hold on;
% for lik = 1:length(unique(alldays(prediction_day_indices).tt(:,3)))
%     for tb = 1:length(firing_perc{1})
%         
%         EFP = nanvar(firing_expect{1}{tb}(like_ind{1}{lik},:),[],2);
%         FP  = nanvar(firing_absolute{1}{tb}(like_ind{1}{lik},:),[],2);
%         
%         DIFFP{1}(lik,tb) = nanmean(FP - EFP);
%         [DIFFP{2}(lik,tb) ,DIFFP{3}(lik,tb)] = boot_bounds(1000,@nanmean,FP-EFP,2.5,97.5);
%     end
%     
%     plot(xsforplot,DIFFP{1}(lik,:),colors_of_plot{lik});
%     patch([xsforplot fliplr(xsforplot)],[DIFFP{2}(lik,:) fliplr(DIFFP{3}(lik,:))],colors_of_plot{lik},...
%         'FaceAlpha',0.25,'EdgeAlpha',0);
% end
