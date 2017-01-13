areas = {'PMd'};%{'PMd','M1'};
half = cell(2,1);

%% Parameters
area = areas{1};

NDIM = 10;
array = eval(sprintf('alldays(1).%s_units',area));
Baseline_block = 1;
Project_blocks = [1 2 3 4];
dsamplesize = 25; % In milliseconds
smoothcutoff = 0.9; % Percentage of spikes that must be filtered
filtsizes = 150;%50:50:500;
filttype = 'causal';
align_col = 17; % 17 for chosen, 13 for correct

t1_t2 = {[6 -200],[10 0]}; % [column time_offset(ms)]

%% Construct smoothed activity profiles for all blocks/neurons
targs = 1:8;
[Actarray,TL] = deal(cell(length(alldays),1));
for task_block = 1:length(alldays)

    tt = alldays(task_block).tt; 

    t_rast = trial_raster(array,tt,t1_t2{1}-[0 3*filtsizes],t1_t2{2}+[0 3*filtsizes]);

    smrast = cell(length(t_rast),1);
    for i = 1:length(t_rast) 
        clc; fprintf('%d - %d/%d\n',task_block, i,length(t_rast));
        for j = 1:size(t_rast{i},1)
            for k = 1:length(filtsizes)
                if strcmp(filttype,'causal')
                    [smthd,mf] = train2cont(t_rast{i}(j,:),filtsizes(k),1);
                else
                    [smthd,mf] = train2cont(t_rast{i}(j,:),filtsizes(k));
                end
                smthd(1:(3*filtsizes)) = []; smthd((end-3*filtsizes+1):end) = [];
                if length(filtsizes) > 1
                    [pv,pl] = findpeaks(smthd);
                    if sum(pv==mf)/sum(t_rast{i}(j,:)) > smoothcutoff
                        break
                    end
                end
            end
            smrast{i}(j,:) = downsample(smthd,dsamplesize);
        end
    end

    % separate by neuron
    smarray = cell(length(array),1);
    for i = 1:length(array)
        smarray{i} = cell2mat(transpose(cellfun(@(x) x(i,:),smrast,'UniformOutput',0)));
    end
    Actarray{task_block} = vertcat(smarray{:})';
    TL{task_block} = cellfun(@(x) size(x,2),smrast);

    clear t_rast
end
clear array
ActarrayO = Actarray; clear Actarray

%% Precompute things that only need to be done once
alcols = { t1_t2{1} ,  [6   0]   ;...
          [6     0] ,  [7   0]   ;...
          [7     0] ,  [8   0]   ;...
          [8     0] ,  [9   0]   ;... 
          [20 -200] , [20   0]   ;...
          [20    0] , t1_t2{2}   };    
t2i = @(t,t1) round(1000*(t-t1)/dsamplesize)+1;
co2t = @(blck,col_off,trial) alldays(blck).tt(trial,col_off(1))+col_off(2)/1000;
INDS = cell(length(alldays),1);
for task_block = 1:length(alldays)
    for i = 1:size(alldays(task_block).tt,1) 
        tstart = alldays(task_block).tt(i,t1_t2{1}(1))+t1_t2{1}(2)/1000; 
        INDS{task_block}(i,:) = reshape(cellfun(@(x) t2i(co2t(task_block,x,i),tstart),alcols)',1,[]);
    end
    INDS{task_block}(INDS{task_block}<0) = NaN;
    INDS{task_block} = mat2cell(INDS{task_block},size(INDS{task_block},1),2*ones(1,length(alcols)));
end
LENS = max(cell2mat(cellfun(@(Y) cellfun(@(x) min(diff(x,[],2)),Y),INDS,'UniformOutput',0)));
CLENS = [1 cumsum(LENS)];
autosplit = cell(length(LENS),1);
for i = 1:length(LENS)
    autosplit{i} = floor(linspace(CLENS(i),CLENS(i+1),floor(LENS(i)./min(LENS))));
    if length(autosplit{i})==1; autosplit{i} = [CLENS(i) autosplit{i}]; end
    autosplit{i} = [autosplit{i}(1:end-1)', autosplit{i}(2:end)',i*ones(length(autosplit{i})-1,1)];
end

ROIs = mat2cell(vertcat(autosplit{:}),ones(sum(cellfun(@(x) size(x,1),autosplit)),1),3);
auto_ROI = cellfun(@(x) x(1:2),ROIs,'UniformOutput',0);

[targids,tinds,cs] = deal(cell(length(alldays),1));
for i = 1:length(Project_blocks)
    bli = Project_blocks(i);
    targids{bli} = mod(round(mod(alldays(bli).tt(:,align_col)+4*pi,2*pi)./(pi/4))+1,9);
    targids{bli}(targids{bli}==0)= 1;
    cs{bli} = [0; cumsum(TL{bli})];
    tinds{bli} = [cs{bli}(1:end-1)+1, cs{bli}(2:end)];
end

TTS = cell(length(Project_blocks),1);
for blck  = 1:length(Project_blocks)
    TTS{blck} = alldays(blck).tt;
end
%%
RS = 500;
Yrep = cell(RS,3);
StaySwap = cell(RS,2);
Lar = size(ActarrayO{1,1},2);
RIS = cell(RS,1);%,Lar);
Biases12 = cell(RS,1);
TAB = cell(RS,2);
%%
parfor reps = 1:RS

    rinds = randperm(Lar);
    RIS(reps) = {rinds};
    for qq = 1:2

        Actarray = cellfun(@(x) x(:,rinds((.5*(qq-1)*Lar + 1):(0.5*qq*Lar))),ActarrayO,'UniformOutput',0);

        %% Do PCA
        [PCax,~,eigenval] = pca(Actarray{Baseline_block});
    %     [PCax,~,eigenval] = pca(vertcat(Actarray{:}));

        %% All activity projected onto 1-target PCs
        traces = cell(length(alldays),1);
        for i = 1:length(Project_blocks)
            bli = Project_blocks(i);
            proj2one = Actarray{bli}*PCax;
            for j = 1:size(tinds{bli},1)
                traces{bli}{j} = proj2one(tinds{bli}(j,1):tinds{bli}(j,2),1:10);
            end
        end

        %% Line up traces for each epoch and store in EPOCH
        [EPOCH] = deal(cell(length(alldays),1));
        for tb = 1:length(Project_blocks)
            bli = Project_blocks(tb); % bli = block
            for i = 1:size(alldays(bli).tt,1) % i = trial
                for j = 1:10 % j = PC
                    for k = 1:size(alcols,1) % k = epoch
                        ep_start = INDS{bli}{k}(i,1);
                        % trim epoch to the minimum size
                        ep_end = min(INDS{bli}{k}(i,2),INDS{bli}{k}(i,1)+LENS(k))-1;

                        if k < size(alcols,1)
                            ep_next = INDS{bli}{k+1}(i,1);
                            ep_end = min(ep_end,ep_next);
                        end

                        if ep_end < ep_start || ep_end > size(traces{bli}{i},1)
                            ep_end = ep_start;
                        end

                        if isnan(ep_start+ep_end)
                           EPOCH{bli}{j,k}(i,:) = nan(1,LENS(k));
                        else
                            padnan = nan(1,LENS(k)-ep_end+ep_start-1);
                            EPOCH{bli}{j,k}(i,:) = [traces{bli}{i}(ep_start:ep_end,j)' padnan];
                        end
                    end

                end
            end
        end

        %% Find N-dimensional directional basis functions 
        DIRTRACES = deal(cell(NDIM,length(INDS)));
        DIRTBT = cell(length(targs),1);
        for j = 1:NDIM % dimension
            for i = 1:length(targs)% target direction
                for k = 1:size(alcols,1) % Epoch
                    DIRTRACES{j,k}(i,:) = nanmean(EPOCH{Baseline_block}{j,k}(targids{Baseline_block}==targs(i),:),1);
                    for t = 1:size(EPOCH{Baseline_block}{j,k},2)
                        DIRTBT{i}{1,k}{1,t}(:,j) = EPOCH{Baseline_block}{j,k}(targids{Baseline_block}==targs(i),t);
                    end
                end
            end
        end
        DIRTBT = cellfun(@(x) horzcat(x{:}),DIRTBT,'UniformOutput',0);
        for i = 1:length(DIRTBT)
            for j = 1:length(DIRTBT{i})
                DIRTBT{i}{j}(isnan(sum(DIRTBT{i}{j},2)),:) = [];
            end
        end
        DIRREG = cell(length(targs),1);
        roi = auto_ROI;

        roicents = cellfun(@mean,roi);
        for i = 1:length(roi)
            for j = 1:length(targs)
                DIRREG{j}{i} = cell2mat(DIRTBT{j}(roi{i})');
            end
        end
        DIRREG = vertcat(DIRREG{:});

        sizecheck = @(x) cellfun(@(x) size(x,1)>size(x,2), x);
        bad_bases = sum(sizecheck(DIRREG)) < size(DIRREG,1);

        DIRREG(:,bad_bases) = [];

        %% Calculate Mahalinobis distances (ALL) 
        bli = Project_blocks(2); % current experiment block
        T = sum(cellfun(@(x) size(x,2),EPOCH{bli}(1,:))); % time points
        P = size(alldays(bli).tt,1); % trials

        % initialize
        
        MD_d = zeros(length(targs),T,P);

        NTRmat = reshape(permute(cell2mat(cellfun(@(x) permute(x',[1 3 2]),EPOCH{bli}','Uni',0)),[1 3 2]),[],10);
        Mahr_rep = cell2mat(cellfun(@(x) reshape(-sqrt(mahal(NTRmat,x)'),1,1,T,P),DIRREG,'Uni',0));
        for i = 1:P
            for t = 1:T
                Maldists = squeeze(Mahr_rep(:,:,t,i));
                maxoverdir = max(Maldists,[],1);
                diratmax = nanmean(Maldists(:,maxoverdir==max(maxoverdir)),2);

                if isempty(diratmax)
                    MD_d(:,t,i) = nan(8,1);
                else
                    MD_d(:,t,i) = circshift(diratmax,5-targids{bli}(i));
                end
            end
        end
        
        if qq == 1
            MD1 = MD_d;
        else
            MD2 = MD_d;
        end

%         for i = 1:P
%             NTR = cell2mat(cellfun(@(x) x(i,:),EPOCH{bli},'Uni',0))';
%             NTRmat = cell2mat(cellfun(@(x) 
%             Mahr = permute(reshape(cell2mat(cellfun(@(x) -sqrt(mahal(NTR,x)'),DIRREG,'UniformOutput',0)),8,T,[]),[1 3 2]);
%             for t = 1:T
%                 Maldists = Mahr(:,:,t);
%                 maxoverdir = max(Maldists,[],1);
%                 diratmax = nanmean(Maldists(:,maxoverdir==max(maxoverdir)),2);
% 
%                 if isempty(diratmax)
%                     MD_d(:,t,i) = nan(8,1);
%                 else
%                     MD_d(:,t,i) = circshift(diratmax,5-targids{bli}(i));
%                 end
% 
%             end
%         end
% 
%         if qq==1
%             MD1 = MD_d;
%         else
%             MD2 = MD_d;
%         end


    end

    mem = 49:63;
    timeall = 1:CLENS(end);
    region = mem;
    trials = [];
    trials.stays{2}.all = find(reshape(nanmean(diff(MD1([1 5],region,:),1),2),[],1)>0);
    trials.swaps{2}.all = find(reshape(nanmean(diff(MD1([1 5],region,:),1),2),[],1)<0);
    for j = 1:8
        trials.stays{2}.targs{j,:} = trials.stays{2}.all(targids{2}(trials.stays{2}.all)==j);
        trials.swaps{2}.targs{j,:} = trials.swaps{2}.all(targids{2}(trials.swaps{2}.all)==j);
    end

    %% Calc diff between targ A and B 
    block = 2;
    trials2use = {'stays','swaps'};% wrongs{4}];%stayswap{2,6};
    target2use = {'all'};%{[1 5],[2 6],[3 7],[4 8]};%{1,2,3,4,5,6,7,8};%{1,2,3,4,5,6,7,8};
    timeper = timeall;

    % trialsets = {corrects{2}(1),wrongs{2}(2)};
    [ys,ub,lb] = deal(cell(length(target2use),length(areas)));
    % [lb, ub] = deal(cell(length(trialsets),1));
    for i = 1:2
        T = [max(timeper(1)-5,0) min(timeper(end)+5,CLENS(end))];
        xs = timeper;
        for k = 1:length(trials2use)
            for j = 1:length(trials2use)
               
                target2use{1} = 1:8;
                trialset = cell2mat(trials.(trials2use{j}){block}.targs(abs(target2use{1})));

                trialsign = cell2mat(cellfun(@(x) repmat(x(1),x(2),1),mat2cell([sign(target2use{1}); ...
                                    cellfun(@(x) size(x,1),trials.(trials2use{j}){block}.targs(abs(target2use{1})))'],2,ones(length(target2use{1}),1)),...
                                    'UniformOutput',0)');

                if i ==1 
                    blockmat = MD1(:,timeper,trialset);
                else
                    blockmat = MD2(:,timeper,trialset);
                end
                blockmat(:,:,trialsign<0) = circshift(blockmat(:,:,trialsign<0),-4,1);
                
                ys{j,i} = nanmean(squeeze(diff(blockmat([1 5],:,:),[],1)),2);
                
                [lb{j,i},ub{j,i}] = boot_bounds(100,@(x) nanmean(x,1),squeeze(diff(blockmat([1 5],:,:),[],1))',2.5,97.5);
            end
            
        end
    end

   
    Yrep(reps,:) = {ys,lb,ub};
    StaySwap(reps,:) = {trials.stays{block}.all, trials.swaps{block}.all};
       
    Biases12(reps) = {[squeeze(nanmean(diff(MD1([1 5],region,:),1),2)), squeeze(nanmean(diff(MD2([1 5],region,:),1),2))]};
    TAB(reps,:) = {squeeze(nanmean(MD1([5 1],region,:)-repmat(nanmean(MD1([3 7],region,:),1),2,1,1),2))', ...
                 squeeze(nanmean(MD2([5 1],region,:)-repmat(nanmean(MD2([3 7],region,:),1),2,1,1),2))'};

    fprintf('repetition: %d\n',reps);
end
%%
[Ymean_stay,Ymean_swap] = deal(cell(2,1));
figure; hold on; 
bnds = [5 95];
for i = 1:2

    Ymean_stay{i} = cell2mat(cellfun(@(x) x{1,i},Yrep(:,1)','UniformOutput',0));
    Ymean_swap{i} = cell2mat(cellfun(@(x) x{2,i},Yrep(:,1)','UniformOutput',0));

    subplot(1,2,i); hold on; 

    plot(nanmean(Ymean_stay{i},2),'b','LineWidth',2); plot(prctile(Ymean_stay{i},bnds,2),'k');
    plot(nanmean(Ymean_swap{i},2),'r','LineWidth',2); plot(prctile(Ymean_swap{i},bnds,2),'k');
    
    Ylimit = ceil(max(abs(reshape(prctile([Ymean_stay{i}; Ymean_swap{i}],bnds,2),[],1))));
    
    plot(repmat(cumsum(LENS),2,1),repmat([-Ylimit;Ylimit],1,length(LENS)),'Color',[.5 .5 .5])
    plot([0 126],[0 0],'k');
end


%%
B12 = cell2mat(Biases12');
B1 = cell2mat(cellfun(@(x) x(:,1),Biases12,'Uni',0)');
B2 = cell2mat(cellfun(@(x) x(:,2),Biases12,'Uni',0)');

[biased,pv] = deal(zeros(size(B12,1),1));
for i = 1:size(B12,1)
    [bv,pv(i)] = ttest(B1(i,:));
    biased(i,:) = bv*sign(mean(B1(i,:)));
end
bias_PNZ = {find(biased>0), find(biased<0), find(biased==0)};

