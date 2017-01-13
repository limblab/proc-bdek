%tune_ranges = { -200:25:1000 , -25:25:500 };
tune_ranges = { [-200:50:1000] , -50:50:450 }; 
%tune_ranges = {[0:50:1000],[0:50:500],[0:50:1000],[0:50:400]};
%tune_ranges = {[0:50:1000],[0:50:1000],[0:50:400]};
%tune_ranges = { [50:50:250 600:50:800] , -50:50:0 }; 
%tune_ranges = {[50 250 600 800],[0 50]};
tune_align = {'target','go'};
%tune_align = {6,7,9};
brain_area = 'PMd';
pred_days = [1 2];

%%
[unit_counts,firing_absolute,spike_count_raw] = deal(cell(length(pred_days),1));
for day=1:length(pred_days)
    prediction_day = pred_days(day);
    for loopthrough = 1:length(tune_align)
        time_bins = tune_ranges{loopthrough};
        time_align = tune_align{loopthrough};

        day_units = eval(sprintf('alldays(1).%s_units',brain_area));
        if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end

        day_pred = prediction_day;
        tt_pred = alldays(day_pred).tt;

        % find rasters for each unit
        unit_pred = cell(length(day_units),1);
        for q = 1:length(day_units)

            clc;
            fprintf('Trial block: %d/%d\nAlignment: %d/%d\nUnit: %d/%d\n',...
                day,length(pred_days),loopthrough,length(tune_align),q,length(day_units));

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
    for loopthrough = 1:length(tune_align)
        time_bins = tune_ranges{loopthrough};
        time_align = tune_align{loopthrough};

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

