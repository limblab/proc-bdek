%-%
brain_area = 'PMd';
Train_block = 1;
Test_block = 2;
 
%-% 
day_units = eval(sprintf('alldays(1).%s_units',brain_area));
if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end

tt_train = alldays(Train_block).tt;
tt_test = alldays(Test_block).tt;

[test_rast] = trial_raster(day_units,tt_test,[6 500],[6 1000]);
[train_rast] = trial_raster(day_units,tt_train,[6 500],[6 1000]);

targ_locs = unique(tt_train(:,13));
find_targ_ind = @(x) find(abs(circ_dist(x,targ_locs))==min(abs(circ_dist(x,targ_locs))));

train_targs = cellfun(find_targ_ind,num2cell(tt_train(:,17)));
test_targs = cellfun(find_targ_ind,num2cell(tt_test(:,17)));
PD_targcell = cellfun(find_targ_ind,num2cell(PDS),'UniformOutput',0);
PD_targs = nan(length(PD_targcell),1); 
for i = 1:length(PD_targs); if ~isempty(PD_targcell{i}); PD_targs(i) = PD_targcell{i}; end; end
OD_targs = mod(PD_targs-1+4,8)+1;
ORTH1_targs = mod(PD_targs-1+2,8)+1;
ORTH2_targs = mod(PD_targs-1-2,8)+1;
[PDmat,ODmat,PDODmat,ORTHmat]= deal(nan(8,length(PD_targs)));
for i = 1:length(PD_targs)
    if ~isnan(PD_targs(i))
        PDmat(PD_targs(i),i) = 1; 
        ODmat(OD_targs(i),i) = 1;
        
        PDODmat(PD_targs(i),i) = 1; PDODmat(OD_targs(i),i)=1;
        
        ORTHmat(ORTH1_targs(i),i) = 1; ORTHmat(ORTH2_targs(i),i) = 1;

    end
end
targ_inds = unique(test_targs);
%%
test_counts = cell2mat(cellfun(@(x) sum(x,2)',test_rast,'UniformOutput',0));
train_counts = cell2mat(cellfun(@(x) sum(x,2)',train_rast,'UniformOutput',0));

%% Loop through target directions and compile training data
targ_means = zeros(length(targ_locs),size(train_counts,2));
for i = 1:length(targ_locs)
    targ_means(i,:) = mean(train_counts(train_targs==i,:));
end
%%
targdist_func = @(x,y) min([mod(x-y,8) mod(y-x,8)]);
%% Loop through testing trials
[pTs,pTs2] = deal(cell(length(test_targs),1));
[logsums,rel_p,logsumsP1,logsumsP2,logsums_aligned] = deal(zeros(length(test_targs),length(targ_locs)));
[logratio,logratioP1,logratioP2,logratioW,sclbest] = deal(zeros(length(test_targs),1));
[t1vt2,rel_ts,rel_tsP1,rel_tsP2,logsumsT] = deal(zeros(length(test_targs),2));
dact = zeros(length(test_targs),size(targ_means,2));
dactP = cell(length(test_targs),3);
nbn_ts = cell(2,1);
pref1 = zeros(length(test_targs),size(targ_means,2));
logsumsPOT = zeros(length(test_targs),3);
for i = 1:length(test_targs)
    targind = test_targs(i);
    oppind = mod(targind+3,8)+1;
    orthind1 = mod(targind+1,8)+1;
    orthind2 = mod(targind-3,8)+1;
    pref1(i,:) = targ_means(targind,:)>targ_means(oppind,:);
    
%     sclbest(i) = fminsearch(@(x) -prob_help_func(test_counts(i,:),x*targ_means,8),1);   
    % Do Bayes calculation
   
    pAgivenTs = poisspdf(repmat(test_counts(i,:),length(targ_locs),1),targ_means);
    pAgivenTs(pAgivenTs==0) = eps;
    
%     pAgivenTs = pAgivenTs.*ORTHmat;
    
    marg = nansum(pAgivenTs,1);
    pTsgivenA = pAgivenTs./repmat(marg,length(targ_locs),1);%.*PDmat;
     
    hfunc = @(x) -prob_help_func(test_counts(i,:),targ_means,8);

    pTs{i} = pTsgivenA;
    pTs2{i} = pTsgivenA([targind oppind],:);
    logsums(i,:) = nansum(log(pTsgivenA),2)';
    logsums_aligned(i,:) = circshift(logsums(i,:),1-test_targs(i),2);

    decoded(i) = find(logsums_aligned(i,:)==max(logsums_aligned(i,:)));
    
    logsumsP1(i,:) = nansum(log(pTsgivenA(:,pref1(i,:)==1)),2)';
    logsumsP2(i,:) = nansum(log(pTsgivenA(:,pref1(i,:)==0)),2)';
    rel_p(i,:) = prod(pTsgivenA,2)'./sum(prod(pTsgivenA,2))';
    logsumsT(i,:) = logsums(i,[targind oppind]);

    logratio(i) = -diff(logsums(i,[targind oppind]));
    logratioP1(i) = -diff(logsumsP1(i,[targind oppind]));
    logratioP2(i) = -diff(logsumsP2(i,[targind oppind]));
%     logratioW(i) = -diff(logsums(i,[targind oppind]).*weightings([targind oppind]));
    
    t1vt2(i,:) = logsums(i,[targind oppind]);
    rel_ts(i,:) = (prod(pTsgivenA([targind oppind],:),2)./sum(prod(pTsgivenA(:,:),2)))';
    rel_tsP1(i,:) = (prod(pTsgivenA([targind oppind],pref1(i,:)==1),2)./sum(prod(pTsgivenA(:,pref1(i,:)==1),2)))';
    rel_tsP2(i,:) = (prod(pTsgivenA([targind oppind],pref1(i,:)==0),2)./sum(prod(pTsgivenA(:,pref1(i,:)==0),2)))';
    nbn_ts{1}(i,:) = pTsgivenA(targind,:); nbn_ts{2}(i,:) = pTsgivenA(oppind,:);
    
%     pdneurs = targ_means(targind,:)==max(targ_means,[],1);
%     odneurs = targ_means(oppind,:)==max(targ_means,[],1);
%     orthneurs = targ_means(orthind1,:)==max(targ_means,[],1) | targ_means(orthind2,:)==max(targ_means,[],1);
%     
    threetargs = mod([targind-1 targind targind+1]-1,8)+1;
    Othreetargs = mod([oppind-1 oppind oppind+1]-1,8)+1;
    
    pdneurs = sum(targ_means(threetargs,:)==repmat(max(targ_means,[],1),3,1));
    odneurs = sum(targ_means(Othreetargs,:)==repmat(max(targ_means,[],1),3,1));
    orthneurs = pdneurs+odneurs == 0;
    
    
    logsumsPOT(i,:) = [-diff(nansum(log(pTsgivenA([targind oppind],pdneurs==1)),2)'),...
                      -diff(nansum(log(pTsgivenA([targind oppind],odneurs==1)),2)'),...
                      -diff(nansum(log(pTsgivenA([targind oppind],orthneurs==1)),2)')];
    
    dact(i,:) = test_counts(i,:)-targ_means(targind,:);
    dactP{i,1} = test_counts(i,pdneurs==1)-targ_means(targind,pdneurs==1);
    dactP{i,2} = test_counts(i,odneurs==1)-targ_means(oppind,odneurs==1);
    dactP{i,3} = dact(i,orthneurs==1);
    
    clc; fprintf('%d/%d\n',i,length(test_targs));
end
