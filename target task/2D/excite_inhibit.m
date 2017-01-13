baselines = nanmean(TA{1},2);

excites = cellfun(@(x) nanmean(x,2)>baselines,TA(2:end),'UniformOutput',0);
excites = [2*ones(size(excites{1})), excites];

inhibits = cellfun(@(x) nanmean(x,2)<baselines,TA(2:end),'UniformOutput',0);
inhibits = [2*ones(size(inhibits{1})), inhibits];