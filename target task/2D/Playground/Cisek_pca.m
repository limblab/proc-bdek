[PCavs, PCall, PCax, Actarray,TL] = deal(cell(3,1));
for task_block = 4:5
    array = alldays(1).PMd_units;
    tt = alldays(task_block).tt; 

    dsamplesize = 25; % Downsample (ms)
    smoothcutoff = 0.10;

    [t_rast] = trial_raster(array,tt,[6 250],[10 0]);
    filtsizes = 50:50:500;

    %%
    smrast = cell(length(t_rast),1);
    for i = 1:length(t_rast) 
        clc; fprintf('%d - %d/%d\n',task_block, i,length(t_rast));
        for j = 1:size(t_rast{i},1)
            for k = 1:length(filtsizes)
                [smthd,mf] = train2cont(t_rast{i}(j,:),filtsizes(k));
                [pv,pl] = findpeaks(smthd);
                if sum(pv==mf)/sum(t_rast{i}(j,:)) < smoothcutoff
                    break
                end
            end
            smrast{i}(j,:) = downsample(smthd,dsamplesize);
        end
    end

    %% separate by neuron
    smarray = cell(length(array),1);
    for i = 1:length(array)
        smarray{i} = cell2mat(transpose(cellfun(@(x) x(i,:),smrast,'UniformOutput',0)));
    end
    totarray = vertcat(smarray{:});
    Actarray{task_block} = totarray';

    [ax, proj] = pca(totarray');

    %% Get inds for each trial
    triallengths = cellfun(@(x) size(x,2),smrast);
    cs = [0; cumsum(triallengths)];
    tinds = [cs(1:end-1)+1, cs(2:end)];

    %% PC traces by trial
    traces = cell(size(tinds,1),1);
    for i = 1:size(tinds,1)
        traces{i} = proj(tinds(i,1):tinds(i,2),1:10);
    end
    %% identify target direction
    targids = round(tt(:,13)./(pi/4))+1;
    targs = unique(targids);
    %% Trim traces to relevant (and equal length) region
    kpinds = 1:(750/dsamplesize+1);
    trimmed = cell(10,1);
    for i = 1:length(traces) % trial
        for j = 1:10 % PC
            trimmed{j}(i,:) = traces{i}(kpinds,j)';
        end
    end
    %% Plot by target direction
    PCdir = cell(10,1);
    for j = 1:10
        for i = 1:length(targs)
            PCdir{j}(i,:) = nanmean(trimmed{j}(targids==targs(i),:),1);
        end
    end
    
    PCavs{task_block} = PCdir;
    PCall{task_block} = trimmed;
    PCax{task_block} = ax;
    
    TL{task_block} = cellfun(@(x) size(x,2),smrast);
    
    clear t_rast
end
%% All activity projected onto 1-target PCs
proj2one = Actarray{task_block}*PCax{1};

%% PC traces by trial
cs = [0; cumsum(TL{task_block})];
tinds = [cs(1:end-1)+1, cs(2:end)];
ONEtraces = cell(size(tinds,1),1);
for i = 1:size(tinds,1)
    ONEtraces{i} = proj2one(tinds(i,1):tinds(i,2),1:10);
end

%% identify target direction
coldir = 17;
targids = round(alldays(3).tt(:,coldir)./(pi/4))+1;
targids1 = round(mod(alldays(1).tt(:,coldir)+2*pi,2*pi)./(pi/4))+1;
targids2 = round(mod(alldays(2).tt(:,coldir)+2*pi,2*pi)./(pi/4))+1;
targids3 = round(mod(alldays(3).tt(:,coldir)+2*pi,2*pi)./(pi/4))+1;
targids4 = round(mod(alldays(4).tt(:,coldir)+2*pi,2*pi)./(pi/4))+1;
targids5 = round(mod(alldays(5).tt(:,coldir)+2*pi,2*pi)./(pi/4))+1;
%% Trim traces to relevant (and equal length) region
kpinds = 1:(750/dsamplesize+1);
ONEtrimmed = cell(10,1);
for i = 1:length(traces) % trial
    for j = 1:10 % PC
        ONEtrimmed{j}(i,:) = ONEtraces{i}(kpinds,j)';
    end
end

BLCK{1} = cellfun(@(x) x(1:size(alldays(1).tt,1),:),ONEtrimmed,'UniformOutput',0);
BLCK{2} = cellfun(@(x) x((size(alldays(1).tt,1)+1):end,:),ONEtrimmed,'UniformOutput',0);
%% Average by target direction
[BLCK1,BLCK2] = deal(cell(10,1));
for j = 1:10
    for i = 1:length(targs)
        BLCK1{j}(i,:) = nanmean(BLCK{1}{j}(targids1==targs(i),:),1);
        BLCK2{j}(i,:) = nanmean(BLCK{2}{j}(targids2==targs(i),:),1);
    end
end

%% Fit
[bd,bdtbt1,bdtbt2,res1,res2,res] = deal(cell(10,1));
for j = 1:10
    for i = 1:length(targs)
        Xcent = [BLCK1{j}(i:end,:) ; BLCK1{j}(1:(i-1),:)]';
%         bd{j}(i,:) = [ones(size(BLCK2{j},2),1), Xcent]\(BLCK2{j}(i,:)');
        bd{j}(i,:) = Xcent\(BLCK2{j}(i,:)');
        res{j}(i,:) = (BLCK2{j}(i,:)' - Xcent*bd{j}(i,:)')';
    end
    
    for i = 1:size(BLCK{2}{1},1)
        Xctbt = [BLCK1{j}(targids2(i):end,:) ; BLCK1{j}(1:(targids2(i)-1),:)]';
%         bdtbt{j}(i,:) = [ones(size(BLCK2{j},2),1), Xctbt]\(BLCK{2}{j}(i,:)');
        bdtbt2{j}(i,:) = Xctbt\(BLCK{2}{j}(i,:)');
        res2{j}(i,:) = (BLCK{2}{j}(i,:)' - Xctbt*bdtbt2{j}(i,:)')';
        
    end
    for i = 1:size(BLCK{1}{1},1)
        Xctbt = [BLCK1{j}(targid1(i):end,:) ; BLCK1{j}(1:(targid1(i)-1),:)]';
%         bdtbt{j}(i,:) = [ones(size(BLCK2{j},2),1), Xctbt]\(BLCK{2}{j}(i,:)');
        bdtbt1{j}(i,:) = Xctbt\(BLCK{1}{j}(i,:)');
        res1{j}(i,:) = (BLCK{1}{j}(i,:)' - Xctbt*bdtbt1{j}(i,:)')';
    end
end

%% Fit all PCs simultaneously
B_cat = zeros(size(BLCK{2}{1},1),8);
for i = 1:size(BLCK{2}{1},1)    
    Xc_cat = cell2mat(cellfun(@(x) [x(targid2(i):end,:) ; x(1:(targid2(i)-1),:)],BLCK1,'UniformOutput',0)')';
    Yc_cat = cell2mat(cellfun(@(x) x(i,:),BLCK{2},'UniformOutput',0)')';
    B_cat(i,:) = Xc_cat\Yc_cat;
end
    
    
%% Identify states associated with each stage and direction
INDS = cell(4,1); % 1:OT   2:MEM   3:CT   4:GO
for i = 1:size(alldays(task_block).tt,1)
    INDS{1}(i,:) = 1:(750/dsamplesize+1);
    INDS{2}(i,:) = round(1000*(diff(alldays(task_block).tt(i,[6 7]))-0.25)/dsamplesize):...
                   (round(1000*(diff(alldays(task_block).tt(i,[6 7]))-0.25)/dsamplesize)+14);
    INDS{3}(i,:) = round(1000*(diff(alldays(task_block).tt(i,[6 8]))-0.25)/dsamplesize):...
                   (round(1000*(diff(alldays(task_block).tt(i,[6 8]))-0.25)/dsamplesize)+10);   % 40
    INDS{4}(i,:) = round(1000*(diff(alldays(task_block).tt(i,[6 9]))-0.25)/dsamplesize):...
                  (round(1000*(diff(alldays(task_block).tt(i,[6 9]))-0.25)/dsamplesize)+12); 
end
% TS = {1:31, 32:46, 47:87 , 88:100};
TS = {1:31, 32:46, 47:57 , 58:70};
%%
EPOCH = cell(10,1);
for i = 1:size(alldays(task_block).tt,1) % trial
    for j = 1:10 % PC

        for k = 1:length(INDS)
            EPOCH{j}{k}(i,:) = ONEtraces{i}(INDS{k}(i,:),j)';
        end
        
    end
end
EPOCH_cat = cellfun(@(x) horzcat(x{:}),EPOCH,'UniformOutput',0);
EPOCHS{1} = cellfun(@(x) x(1:size(alldays(1).tt,1),:),EPOCH_cat,'UniformOutput',0);
EPOCHS{2} = cellfun(@(x) x((size(alldays(1).tt,1)+1):end,:),EPOCH_cat,'UniformOutput',0);

%%
[EPOCH1,EPOCH2] = deal(cell(10,1));
for j = 1:10
    for i = 1:length(targs)
        EPOCH1{j}(i,:) = nanmean(EPOCHS{1}{j}(targids1==targs(i),:),1);
        EPOCH2{j}(i,:) = nanmean(EPOCHS{2}{j}(targids2==targs(i),:),1);
    end
end
BYDIR = cell(8,1);
for i = 1:length(BYDIR)
    BYDIR{i} = cell2mat(cellfun(@(x) x(i,:),EPOCH1,'UniformOutput',0));
end

BYTRI = cell(size(EPOCHS{2}{1},1),1);
for i = 1:size(EPOCHS{2}{1},1)
    BYTRI{i} = cell2mat(cellfun(@(x) x(i,:),EPOCHS{2},'UniformOutput',0));
end

EPOCH_BL = EPOCH1;

%% Fit all PCs per time point simultaneously
B_time = zeros(8,size(EPOCHS{2}{1},2),size(EPOCHS{2}{1},1));
for i = 1:size(EPOCHS{2}{1},1)    
    for t = 1:size(EPOCHS{2}{1},2)
        Xc_tim = cell2mat(cellfun(@(x) [x(targids5(i):end,t) ; x(1:(targids5(i)-1),t)],EPOCH_BL,'UniformOutput',0)')';
        Yc_tim = cell2mat(cellfun(@(x) x(i,t),EPOCHS{2},'UniformOutput',0)')';
        
        B_time(:,t,i) = Xc_tim\Yc_tim;
    end
end
Bt_av = mean(B_time,3);

%% Fit each epoch separately
B_epoch = zeros(8,length(TS),size(EPOCHS{2}{1},1));
for i = 1:size(EPOCHS{2}{1},1) % Trials
    for ep = 1:length(TS) % Epoch
        Xc_tim = cell2mat(cellfun(@(x) [x(targids5(i):end,TS{ep}) ; x(1:(targids5(i)-1),TS{ep})],EPOCH_BL,'UniformOutput',0)')';
        Yc_tim = cell2mat(cellfun(@(x) x(i,TS{ep}),EPOCHS{2},'UniformOutput',0)')';
        
        B_epoch(:,ep,i) = Xc_tim\Yc_tim;
    end
end
Be_av = mean(B_epoch,3);
%%
Be_avD = cell(length(targs),1);
for i = 1:length(targs)
    Be_avD{i} = mean(B_epoch(:,:,targids2==targs(i)),3);
end

%%
for i = 1:length(BYTRI)
    p = BYTRI{i};
    [dst,phse] = deal(zeros(8,size(ds,1)));
    for j = 1:8
        ds = pdist2(p',BYDIR{j}');
        ws = 1-(ds-min(ds(:)))./(max(ds(:))-min(ds(:)));
        for k = 1:size(ds,1)
            [dst(j,k), phse(j,k)] = min(ws(k,:));
        end
    end
end
    