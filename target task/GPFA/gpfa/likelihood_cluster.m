%% concatenated trajs for each likelihood condition
[us,liketag] = deal(cell(length(seqlikes),1));
for i = 1:length(seqlikes)  
    us{i} = [];
    for j = 1:length(seqlikes{i})  
        us{i} = [us{i}; seqlikes{i}(j).xorth(1:3,:)'];
    end
    liketag{i} = i*ones(length(us{i}),1);
end
us_all = vertcat(us{:});
tag_all = vertcat(us{:});

%%
dim = [2 3];

idx = kmeans(us_all(:,dim),length(us));

unlengths = cellfun(@length,us);

[idx_like,us_tri,trilengthsk,idx_tri,trilengths,corr_ts,corr_tri,...
    classav,correct_class,class_times,corr_times] = deal(cell(length(seqlikes),1));

idx_like{1} = idx(1:unlengths(1),:);
for i = 2:length(seqlikes)
    idx_like{i} = idx((unlengths(i-1)+1):sum(unlengths(1:i)),:);
end

for i = 1:length(seqlikes)

    trilengths{i} = horzcat(seqlikes{i}.T);
    
    us_tri{i}{1} = us{i}(1:trilengths{i}(1),:);
    idx_tri{i}{1} = idx_like{i}(1:trilengths{i}(1));
    corr_times{i}{1} = idx_tri{i}{1}==i;
    corr_ts{i}{1} = find(corr_times{i}{1});
    corr_tri{i}{1} = length(corr_ts{i}{1})./length(idx_tri{i}{1});

    for j = 2:length(trilengths{i})
        
        trial_inds = (sum(trilengths{i}(1:(j-1)))+1):sum(trilengths{i}(1:j));
        
        us_tri{i}{j} = us{i}(trial_inds,:);
        idx_tri{i}{j} = idx_like{i}(trial_inds,:);
        
        corr_times{i}{j} = idx_tri{i}{j}==i;
        corr_ts{i}{j} = find(corr_times{i}{j});
        corr_tri{i}{j} = length(corr_ts{i}{j})./length(idx_tri{i}{j});
    end
    
    classav{i} = vertcat(corr_tri{i}{:});
    correct_class{i} = classav{i} > 0.5;
    
    class_times{i} = vertcat(corr_ts{i}{:});
end


        