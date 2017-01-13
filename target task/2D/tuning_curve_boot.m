tune_ranges = [0 400 800];
tune_align = 'target';
brain_area = 'PMd';

neurboot = cell(length(tune_ranges)-1,1);
num_boots = 10;
for ranges = 1:length(tune_ranges)-1
    tune_range = tune_ranges(ranges:ranges+1);
    % Initialize
    co_units = eval(sprintf('alldays(1).%s_units',brain_area));
    if strcmp(brain_area,'M1'); co_units(end-1:end) =[]; end
    %neurboot = cell(length(co_units),length(tune_range)-1);

    for i = 1:length(co_units) % Loop through units
         clc; fprintf('%d/%d Tuning: (%d/%d)\n',ranges,length(tune_ranges)-1,i,length(co_units));

        % set time bin edges
        t1 = tune_range(1);
        t2 = tune_range(2);

        [cents,tune] = co_tuning_boot(co_units,alldays(1).tt(:,10),alldays(1).tt,t1,t2,tune_align,i,num_boots);

        % Smooth and fix edges on tuning curves

        neurboot{ranges}{i}.tuning = nan(num_boots,628);
        for z = 1:num_boots
            tune_cents = [cents-2*pi cents cents+2*pi];
            tune_pad = repmat(tune(z,:)',3,1);

            interp_cents = tune_cents(1):.01:tune_cents(end);
            tune_interp = interp1(tune_cents,tune_pad,interp_cents);

            smooth_tune = smooth(tune_interp,100);

            wrapped_cents = interp_cents(interp_cents >=0 & interp_cents< 2*pi);
            wrapped_tune = smooth_tune(interp_cents >= 0 & interp_cents <2*pi); 

            neurboot{ranges}{i}.tuning(z,:) = wrapped_tune;
        end

    end
end
centers = round(wrapped_cents*1000)./1000;

%% obtain spike counts

time_bins = tune_ranges;
time_align = tune_align;
prediction_day_indices = [2];
prior_mean = 90;

xsforplot = .5*(time_bins(1:end-1)+time_bins(2:end));
 
unit_counts = cell(length(prediction_day_indices),1);
for z = 1:length(prediction_day_indices)

    prediction_day = prediction_day_indices(z);
    day_units = eval(sprintf('alldays(%d).%s_units',prediction_day,brain_area));
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

%% Get expected and actual trial counts
trial_counts = cell(length(prediction_day_indices),1);
for bin = 1:length(time_bins) - 1
    clc; fprintf('bin %d/%d\n',bin,length(time_bins)-1);
    for z = 1:length(prediction_day_indices)
        day_units = eval(sprintf('alldays(1).%s_units',brain_area));
        if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end
        for i = 1:size(alldays(prediction_day_indices(z)).tt,1)    
            for q = 1:length(day_units)       
                trial_counts{z}{bin}(i,q) = unit_counts{z}{q}(i,bin);
            end
        end
    end
end
