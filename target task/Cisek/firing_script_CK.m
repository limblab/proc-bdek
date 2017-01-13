%tune_ranges = { -200:25:1000 , -25:25:500 };
%tune_ranges = { -200:50:1500, -50:50:500, -50:50:1500, -50:50:500 ,0:50:500}; 
%tune_ranges = { -200:50:1500, -50:50:500, -50:50:1500, -50:50:500 }; 
%tune_ranges = { -200:50:1000, 0:50:1000 , 0:50:500 }; 
tune_ranges = { -200:50:1000 , -50:50:300, -50:50:1000};
%tune_ranges = { -400:200:1600, -200:200:600, -200:200:600};
%tune_ranges = {[250 800],[0 250]};
tune_align = {6,7,8};
%tune_align = {6, 7, 9};
%tune_align = {6,7,8,9,10};
brain_area = 'PMd';
session_segs = {1,2};

%%
[unit_counts,firing_absolute,spike_count_raw] = deal(cell(length(session_segs),1));
for sess = 1:length(session_segs)
    
    pred_days = session_segs{sess};
    if sess == 1
        tt_pred = vertcat(alldays(pred_days).tt);
    else
        tt_pred = vertcat(alldays(pred_days).tt);
    end
    
    for loopthrough = 1:length(tune_align)
        time_bins = tune_ranges{loopthrough};
        time_align = tune_align{loopthrough};

        day_units = eval(sprintf('alldays(1).%s_units',brain_area));
        if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end

        % find rasters for each unit
        unit_pred = cell(length(day_units),1);
        for q = 1:length(day_units)

            clc;
            fprintf('Loop: %d/%d\nUnit: %d/%d\n',loopthrough,length(tune_align),q,length(day_units));

            [rast_out,rast_inds] = raster_plotCK(day_units{q},tt_pred,[time_bins(1)./1000 time_bins(end)./1000],...
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

        [trial_counts] = deal(cell(length(pred_days),1));
        for bin = 1:length(time_bins) - 1
            clc; fprintf('Loop: %d/%d\nbin %d/%d\n',loopthrough,length(tune_align),bin,length(time_bins)-1);

            day_units = eval(sprintf('alldays(1).%s_units',brain_area));
            if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end

            for i = 1:size(tt_pred,1)   
                for q = 1:length(day_units)       
                    trial_counts{bin}(i,q) = unit_counts{loopthrough}{q}(i,bin);
                end
            end
            firing_absolute{sess}{loopthrough}{bin} = trial_counts{bin};
            spike_count_raw{sess}{loopthrough}{bin} = trial_counts{bin}*(diff(time_bins(1:2))/1000);
        end

    end
end
