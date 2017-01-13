%-%
brain_area = 'PMd';
Train_block = 1;
Test_block = 2;

timewins = 0:50:1000;
alignment = 6*ones(1,length(timewins));

%-% 
day_units = eval(sprintf('alldays(1).%s_units',brain_area));
if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end

tt_train = alldays(Train_block).tt;
tt_test = alldays(Test_block).tt;

targ_locs = unique(tt_train(:,13));
find_targ_ind = @(x) find(abs(circ_dist(x,targ_locs))==min(abs(circ_dist(x,targ_locs))));

train_targs = cellfun(find_targ_ind,num2cell(tt_train(:,17)));
test_targs = cellfun(find_targ_ind,num2cell(tt_test(:,17)));
opp_targs = cellfun(find_targ_ind,num2cell(tt_test(:,17)+pi));
targ_inds = unique(test_targs);

[test_rast,train_rast,test_counts,train_counts,targ_means,logratio,dact,...
    dactPEAK,dactVAL,logsums,logsums12,T12] =...
    deal(cell((length(timewins)-1),1));
for i = 1:(length(timewins)-1)
    [test_rast{i}] = trial_raster(day_units,tt_test,[alignment(i) timewins(i)],[alignment(i) timewins(i+1)]);
    [train_rast{i}] = trial_raster(day_units,tt_train,[alignment(i) timewins(i)],[alignment(i) timewins(i+1)]);
    
    test_counts{i} = cell2mat(cellfun(@(x) sum(x,2)',test_rast{i},'UniformOutput',0));
    train_counts{i} = cell2mat(cellfun(@(x) sum(x,2)',train_rast{i},'UniformOutput',0));
    
    for j = 1:length(targ_locs)
        targ_means{i}(j,:) = mean(train_counts{i}(train_targs==j,:));
    end
    
    
    for k = 1:length(test_targs)
        targind = test_targs(k);
        oppind = mod(targind+3,8)+1;

        pAgivenTs = poisspdf(repmat(test_counts{i}(k,:),length(targ_locs),1),targ_means{i});
        pAgivenTs(pAgivenTs==0) = eps;
        marg = nansum(pAgivenTs,1);
        pTsgivenA = pAgivenTs./repmat(marg,length(targ_locs),1);
        
        pAgivenT12 = pAgivenTs([targind oppind],:);
        PT12givenA = pAgivenT12./repmat(nansum(pAgivenT12,1),2,1);
        
        T12{i}{1}(k,:) = PT12givenA(1,:);
        T12{i}{2}(k,:) = PT12givenA(2,:);

        logsums{i}(k,:) = nansum(log(pTsgivenA),2)';
        logsums12{i}(k,:) = logsums{i}(k,[targind oppind]);
        logratio{i}(k,:) = -diff(logsums{i}(k,[targind oppind]));
        
        dact{i}(k,:) = test_counts{i}(k,:)-targ_means{i}(targind,:);
        dactPEAK{i}(k,:) = test_counts{i}(k,:)-max(targ_means{i},[],1);
        dactVAL{i}(k,:) = test_counts{i}(k,:)-min(targ_means{i},[],1);

        clc; fprintf('%d/%d\n%d/%d\n',i,length(timewins)-1,k,length(test_targs));
    end
    
end

%%
for i = 1:length(T12)
    
    Targ1_for = T12{i}{1}>0.5;
    Targfor{1,1}(:,i) = nansum(T12{i}{1}.*Targ1_for,2)./sum(Targ1_for,2);
     
    Targ2_for = T12{i}{2}>0.5;
    Targfor{2,1}(:,i) = nansum(T12{i}{2}.*Targ2_for,2)./sum(Targ2_for,2);
    
    Targ1_against = T12{i}{1}<0.5;
    Targfor{1,2}(:,i) = nansum(T12{i}{1}.*Targ1_against,2)./sum(Targ1_against,2);
    
    Targ2_against = T12{i}{2}<0.5;
    Targfor{2,2}(:,i) = nansum(T12{i}{2}.*Targ2_against,2)./sum(Targ2_against,2);
end
    
%%
d2rch = @(PD,rch) abs(circ_dist(PD,rch));
[actP] = cell(2,2);
actM = cell(2,1);
for trl = 1:size(test_counts{1},1)
    
    neursPD = d2rch(PDS,alldays(2).tt(trl,17));
    neursPD(neursPD > pi/4) = nan; neursPD(neursPD <= pi/4) = 1;
    
    neursOD = d2rch(PDS,alldays(2).tt(trl,17)+pi);
    neursOD(neursOD > pi/4) = nan; neursOD(neursOD <= pi/4) = 1;
    
    ttarg = test_targs(trl);
    otarg = opp_targs(trl);
    
    for time = 1:length(test_counts)
        
        actP{1,1}(trl,time) = nanmean(test_counts{time}(trl,neursPD==1) - targ_means{time}(ttarg,neursPD==1));
        actP{1,2}(trl,time) = nanmean(test_counts{time}(trl,neursPD==1) - targ_means{time}(otarg,neursPD==1));
        actP{2,1}(trl,time) = nanmean(test_counts{time}(trl,neursOD==1) - targ_means{time}(ttarg,neursOD==1));
        actP{2,2}(trl,time) = nanmean(test_counts{time}(trl,neursOD==1) - targ_means{time}(otarg,neursOD==1));
        
%         actM{1}(trl,time) = nanmean(test_counts{time}(trl,neursPD==1) - mean(targ_means{time}(:,neursPD==1)));
%         actM{2}(trl,time) = nanmean(test_counts{time}(trl,neursOD==1) - mean(targ_means{time}(:,neursOD==1)));
       
%         modlds = max(targ_means{time})-min(targ_means{time});
%         actM{1}(trl,time) = nanmean((test_counts{time}(trl,neursPD==1) - mean(targ_means{time}(:,neursPD==1)))./modlds(neursPD==1));
%         actM{2}(trl,time) = nanmean((test_counts{time}(trl,neursOD==1) - mean(targ_means{time}(:,neursOD==1)))./modlds(neursOD==1));
        
        actM{1}(trl,time) = nanmean((test_counts{time}(trl,neursPD==1) - mean(targ_means{time}(:,neursPD==1)))>0);
        actM{2}(trl,time) = nanmean((test_counts{time}(trl,neursOD==1) - mean(targ_means{time}(:,neursOD==1)))>0);
    end
end

    