tune_range = [500 750];
tune_align = 'target';
brain_area = 'PMd';

% Initialize
co_units = eval(sprintf('alldays(1).%s_units',brain_area));
if strcmp(brain_area,'M1'); co_units(end-1:end) =[]; end
neurons = cell(length(co_units),length(tune_range)-1);

t1 = tune_range(1);
t2 = tune_range(2);
spike_counts = cell(length(alldays),1);
for j = 1:length(alldays) %DAYS
    for i = 1:length(co_units) % Loop through units
        
        cur_tt = alldays(j).tt;
        set_of_likes = unique(cur_tt(:,3));
        
        for zz = 1:length(set_of_likes)
            inds{j}{zz} = find(cur_tt(:,3)==set_of_likes(zz));
        end
    
        indices = inds{j};
        
        for q = 1:length(indices) % likelihood uncertainty
            clc; fprintf('Counting: \nDay %d\nUnit %d/%d\nEpoch %d\n',j,i,length(co_units),q);

            % set time bin edges
            t1 = tune_range(1);
            t2 = tune_range(2);

            [poss,counts] = co_tuning(alldays(j).PMd_units,alldays(j).tt(indices{q},10),alldays(j).tt(indices{q},:),t1,t2,tune_align,i);

            % Smooth and fix edges on tuning curves

            spike_counts{j}{q}(i,:) = counts;

            %{session portion}{trial portion}{unit}
            %neurons{i}.raster = tune_rast;

        end
    end
end
%%
for i = 1:length(spike_counts)-1
    modulations = repmat((max(spike_counts{1}{1},[],2)-min(spike_counts{1}{1},[],2)),1,8);
    good_neurons = find(modulations(:,1) > 0);
   
    for j = 1:length(inds{i+1})
        scd{i}{j} = spike_counts{i+1}{j}(good_neurons,:) - spike_counts{1}{1}(good_neurons,:);
        scd{i}{j} = scd{i}{j}./modulations(good_neurons,:);
        scd{i}{j}(isnan(scd{i}{j})) = 0;
        av_diff{i}{j} = nanmean(scd{i}{j},2);
    end
    
end

%%
for j = 1:length(inds{2})

    [n{j},b{j}] = hist(av_diff{1}{j},25);
    
end
colors_for_plot = {'r','g','b','c'};
figure; hold on;
for j = 1:length(inds{2})
    subplot(1,length(inds{2}),length(inds{2})-j+1);
    bar(b{j},n{j},colors_for_plot{j}); %xlim([-10 10]);
end

%%
for j = 1:length(inds{2})
    bb{j} = skewness(av_diff{1}{length(inds{2})-j+1});
    [lb{j}, ub{j}] = boot_bounds(10000,@skewness,av_diff{1}{length(inds{2})-j+1},2.5,97.5);
end
figure; errorbar(1:length(inds{2}),cell2mat(bb),cell2mat(bb)-cell2mat(lb),cell2mat(ub)-cell2mat(bb));

