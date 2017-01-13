% tune_ranges = [-200:25:800];
tune_ranges = [550:50:850];
tune_align = 'target';
brain_area = 'M1';
prediction_day_indices = 2;

%%
time_bins = tune_ranges;
time_align = tune_align;

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
    end
end
%%
[trial_counts,firing_absolute,spike_count_raw] = deal(cell(length(prediction_day_indices),1));
for bin = 1:length(time_bins) - 1
    clc; fprintf('bin %d/%d\n',bin,length(time_bins)-1);
    for z = 1:length(prediction_day_indices)
        
        day_pred = prediction_day_indices(z);
        tt_pred = alldays(day_pred).tt;  

        day_units = eval(sprintf('alldays(1).%s_units',brain_area));
        if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end
 
        for i = 1:size(tt_pred,1)   
            for q = 1:length(day_units)       
                trial_counts{z}{bin}(i,q) = unit_counts{z}{q}(i,bin);
            end
        end
        firing_absolute{z}{bin} = trial_counts{z}{bin};
        spike_count_raw{z}{bin} = trial_counts{z}{bin}*(diff(time_bins(1:2))/1000);
    end
end

%xsforplot = .5*(time_bins(1:end-1)+time_bins(2:end));
xsforplot = 1;
