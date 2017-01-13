%-%
brain_area = 'PMd';
Train_block = 7;
Test_block = 9;
 
%-% 
day_units = eval(sprintf('alldays(1).%s_units',brain_area));
if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end

tt_train = alldays(Train_block).tt;
tt_test = alldays(Test_block).tt;
tt_test(:,20) = 2; tt_test(tt_test(:,14) > 500,20) = 1;

[test_rast] = trial_raster(day_units,tt_test,[8 -1000],[8 1000]);
[train_rast] = trial_raster(day_units,tt_train,[6 500],[6 1000]);

targ_locs = unique(tt_train(:,13));
find_targ_ind = @(x) find(abs(circ_dist(x,targ_locs))==min(abs(circ_dist(x,targ_locs))));
sliding_window = @(x) conv(x,ones(1,250),'same');
averaging_window = @(x) sum(x)/2;

train_targs = cellfun(find_targ_ind,num2cell(tt_train(:,17)));
test_targs = cellfun(find_targ_ind,num2cell(tt_test(:,17)));
%%
average = zeros(length(day_units),size(tt_train,1));
for i = 1:size(tt_train,1)
    average(:,i) = cellfun(averaging_window,num2cell(train_rast{i},2));
end
%%
[trial_ratios, Irat] = deal(cell(size(tt_test,1),1));
PDrank = nan(size(tt_test,1),length(day_units));
xopt = zeros(size(tt_test,1),2);
for i = 1:size(tt_test,1)
    
    slide = cell2mat(cellfun(sliding_window,num2cell(test_rast{i},2),'UniformOutput',0));
    %slide = sum(test_rast{i},2);
        
    train_trials_for = find(train_targs == find_targ_ind(tt_test(i,17)));
    train_trials_against = find(train_targs == find_targ_ind(tt_test(i,17)+pi));
    
    foraverage = mean(average(:,train_trials_for),2);
    againstaverage = mean(average(:,train_trials_against),2);
    
    valids = find(prod([foraverage againstaverage],2)~=0);
    
    probs_for = poisspdf(slide(valids,:),repmat(foraverage(valids,:),1,size(slide,2)));
    probs_against = poisspdf(slide(valids,:),repmat(againstaverage(valids,:),1,size(slide,2)));

    log_ratios = log(probs_for./probs_against);
    
%     pdval = PD1{1}(valids);
%     dfpd = abs(dot(repmat([cos(tt_test(i,17)) sin(tt_test(i,17))],length(valids),1),[cos(pdval) sin(pdval)],2));
%     
%     [~,utilrank] = sortrows(log_ratios);
%     PDrank(i,1:length(utilrank)) = valids(utilrank);
    
    Irat{i} = log_ratios;
    trial_ratios{i} = sum(log_ratios,1);
    
    
%     helper = @(x) closest_targ(slide,foraverage,againstaverage,x);
%     xopt(i,:) = fminsearch(helper,[1 1]);
%     
%     rat_opt(i) = sum(log((probs_for.*xopt(i,1)+xopt(i,2))./(probs_against.*xopt(i,1)+xopt(i,2))));
    clc; fprintf('Trial: %d/%d\n',i,size(tt_test,1));
end
%%
trial_lengths = cellfun(@length,trial_ratios);
trial_matrix = nan(size(tt_test,1),max(trial_lengths));
for i = 1:size(tt_test,1)
    trial_matrix(i,1:trial_lengths(i)) = trial_ratios{i};
    %trial_matrix(i,(end-trial_lengths(i)+1):end) = trial_ratios{i};
end
mincol = min(min(trial_matrix));
maxcol = max(max(trial_matrix));

cueinds = round(1000*(tt_test(:,9)-tt_test(:,6)));

blocks = unique(tt_test(:,20));
[tlengths,tratios,ord,inds] = deal(cell(length(blocks),1));
figure; hold on; 
for i = 1:length(blocks)
    inds{i} = find(tt_test(:,20)==blocks(i));
    tlengths{i} = trial_lengths(inds{i});
    tratios{i} = trial_matrix(inds{i},:);
    
    %[~,ord{i}] = sortrows(tlengths{i});
    %[~,ord{i}] = sortrows(max(tratios{i},[],2));
    [~,ord{i}] = sortrows(max(tratios{i}(:,500:1500),[],2));
    subplot(1,length(blocks),i); imagesc(tratios{i}(ord{i},:),[mincol,maxcol]); hold on; 
    title(sprintf('Block: %d',blocks(i)),'FontSize',15);
    
    
%     for j = 1:length(inds)
%         plot(cueinds(inds(ord{i}(j))),j,'k.','MarkerSize',20);
%         plot(cueinds(inds(ord{i}(j)))+RT(inds(ord{i}(j))),j,'w.','MarkerSize',20);
%     end
    
end
    