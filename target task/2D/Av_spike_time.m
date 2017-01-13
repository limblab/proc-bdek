tune_ranges = [0 800];
tune_align = 'target';
brain_area = 'PMd';

neuronsT = cell(length(tune_ranges)-1,1);
CO_index = 1;
for ranges = 1:length(tune_ranges)-1
    tune_range = tune_ranges(ranges:ranges+1);
    % Initialize
    co_units = eval(sprintf('alldays(%d).%s_units',CO_index,brain_area));
    if strcmp(brain_area,'M1'); co_units(end-1:end) =[]; end
    %neurons = cell(length(co_units),length(tune_range)-1);

    for i = 1:length(co_units) % Loop through units
         clc; fprintf('%d/%d Tuning: (%d/%d)\n',ranges,length(tune_ranges)-1,i,length(co_units));

        % set time bin edges
        t1 = tune_range(1);
        t2 = tune_range(2);

        [cents,~,~,~,~,AVTIME] = co_tuning(co_units,alldays(CO_index).tt(:,10),...
            alldays(CO_index).tt,t1,t2,tune_align,i);

        % Smooth and fix edges on tuning curves
        front_cent = cents;
        back_cent = cents;

        tune = AVTIME;
        front_tune = AVTIME; 
        back_tune = AVTIME;

        tune_cents = [back_cent-2*pi cents front_cent+2*pi];
        tune_pad = [back_tune;tune;front_tune];

        interp_cents = tune_cents(1):.01:tune_cents(end);
        tune_interp = interp1(tune_cents,tune_pad,interp_cents);

        smooth_tune = smooth(tune_interp,100);

        wrapped_cents = interp_cents(interp_cents >=0 & interp_cents< 2*pi);
        wrapped_tune = smooth_tune(interp_cents >= 0 & interp_cents <2*pi);
        
        neuronsT{ranges}{i}.tuning = wrapped_tune;
    end
end

%%
time_bins = tune_ranges;
time_align = tune_align;
prediction_day_indices = 2;

xsforplot = .5*(time_bins(1:end-1)+time_bins(2:end));
 
[unit_times, like_ind, avtimes, unit_deltas, avdelta] = deal(cell(length(prediction_day_indices),1));

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

        [rast_out,rast_inds,~,AvTime] = raster_plot(day_units{q},tt_pred,[time_bins(1)./1000 time_bins(end)./1000],...
            time_align,0,'none');

        out = zeros(size(vertcat(AvTime{:})));
        for m = 1:length(rast_out)
            out(rast_inds{m},:) = AvTime{m};
        end
        out(out==0)=NaN;

        unit_times{z}{q} = out;
        
        movelocs = round(100*tt_pred(:,10));
        movelocs(movelocs==0) = 1;
        
        eout = neuronsT{1}{q}.tuning(movelocs);
        
        delta_time = out - eout;
        
        unit_deltas{z}{q} = out - eout; 

        set_of_inds = flipud(unique(tt_pred(:,3)));
        for qq = 1:length(set_of_inds)
            like_ind{z}{qq} = find(tt_pred(:,3)==set_of_inds(qq));
            avtimes{z}(q,qq) = nanmean(out(like_ind{z}{qq}));
            avdelta{z}(q,qq) = nanmean(delta_time(like_ind{z}{qq}));
        end
    end
end
%%
% figure; hold on;
% avtimes{1}(avtimes{1} == 0) = NaN;
% TIME = cell(3,1);
% cop = {'b','g','r'};
% for i = 1:size(avtimes{1},2)
%     
%     [TIME{i}(1), TIME{i}(2)] = boot_bounds(1000,@nanmean,avtimes{1}(:,i),2.5,97.5);
%     TIME{i}(3) = nanmean(avtimes{1}(:,i));
%     
% %     plot(i,TIME{i}(3),'.','Color',cop{i},'MarkerSize',15);
% %     plot([i i],TIME{i}(1:2),cop{i},'LineWidth',3);
%     
%     plot(TIME{i}(3),i,'.','Color',cop{i},'MarkerSize',15);
%     plot(TIME{i}(1:2),[i i],cop{i},'LineWidth',3);
%     xlim([-200 1000]); ylim([0 50]);
% end
% title('Average Spike Time (s)');
% 
% %%
% figure; hold on;
% avtimes{1}(avtimes{1} == 0) = NaN;
% DTIME = cell(3,1);
% cop = {'b','g','r'};
% for i = 1:size(avtimes{1},2)
%     
%     [DTIME{i}(1), DTIME{i}(2)] = boot_bounds(1000,@nanmean,avdelta{1}(:,i),2.5,97.5);
%     DTIME{i}(3) = nanmean(avdelta{1}(:,i));
%     
%     plot(i,DTIME{i}(3),'.','Color',cop{i},'MarkerSize',15);
%     plot([i i],DTIME{i}(1:2),cop{i},'LineWidth',3);
% end
% title('Change in Average Spike Time (s)');

%%
timing_array = horzcat(unit_times{1}{:});
delta_array = horzcat(unit_deltas{1}{:});

set_of_inds = flipud(unique(tt_pred(:,3)));
[av_alltimes, l_alltimes, u_alltimes, T_change, ...
    av_alldeltas, l_alldeltas, u_alldeltas] = deal(cell(length(set_of_inds),1));

for qq = 1:length(set_of_inds)
    like_ind{z}{qq} = find(tt_pred(:,3)==set_of_inds(qq));
    
    %%%% Timing %%%%
    alltimes = vertcat(timing_array(like_ind{z}{qq},:));
    alltimes = alltimes(:); alltimes(isnan(alltimes)) = [];
    
    av_alltimes{qq} = mean(alltimes);
    [l_alltimes{qq}, u_alltimes{qq}] = boot_bounds(1000,@mean,alltimes,2.5,97.5);

    %%%% Changes %%%%
    alldeltas = vertcat(delta_array(like_ind{z}{qq},:));
    alldeltas = alldeltas(:); alldeltas(isnan(alldeltas)) = [];
    
    av_alldeltas{qq} = mean(alldeltas);
    [l_alldeltas{qq},u_alldeltas{qq}] = boot_bounds(1000,@mean,alldeltas,2.5,97.5);
    
    T_change{qq} = alldeltas;

end

%%
% figure; hold on;
% subplot(1,2,1);
% for qq = 1:length(set_of_inds)
%     like_ind{z}{qq} = find(tt_pred(:,3)==set_of_inds(qq));
% 
%     subplot(1,2,1); hold on;
%     plot(av_alltimes{qq},qq,'.','Color',cop{qq},'MarkerSize',15);
%     plot([l_alltimes{qq} u_alltimes{qq}],[qq qq],cop{qq},'LineWidth',3);
%     
% end
% subplot(1,2,2); hold on;
% for qq = 1:length(set_of_inds)
%     like_ind{z}{qq} = find(tt_pred(:,3)==set_of_inds(qq));
% 
%     plot(av_alldeltas{qq},qq,'.','Color',cop{qq},'MarkerSize',15);
%     plot([l_alldeltas{qq} u_alldeltas{qq}],[qq qq],cop{qq},'LineWidth',3);
%     
% end
