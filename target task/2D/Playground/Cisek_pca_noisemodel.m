% find_RT; % get reaction times for each block
% for task_block = 1:length(alldays)
%     alldays(task_block).tt(:,20) = alldays(task_block).RT; %#ok<SAGROW>
%     badtrls = find(isnan(alldays(task_block).tt(:,20)));
%     alldays(task_block).tt(badtrls,20) = alldays(task_block).tt(badtrls,9);
%         
% end
areas = {'PMd','M1'};%{'PMd','M1'};

for ba_loop = 1:length(areas)
    %% Parameters
    area = areas{ba_loop};

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

    %% Do PCA
    [PCax,~,eigenval] = pca(Actarray{Baseline_block});
%     [lambda, psi] = factoran(Actarray{Baseline_block},10);
%     [PCax,~,eigenval] = pca(vertcat(Actarray{:}));

    %% All activity projected onto 1-target PCs
    [targids,traces] = deal(cell(length(alldays),1));
    for i = 1:length(Project_blocks)

        bli = Project_blocks(i);

        targids{bli} = mod(round(mod(alldays(bli).tt(:,align_col)+4*pi,2*pi)./(pi/4))+1,9);
        targids{bli}(targids{bli}==0)= 1;
        proj2one = Actarray{bli}*PCax;
%         proj2one = Actarray{bli}*lambda;

        cs = [0; cumsum(TL{bli})];
        tinds = [cs(1:end-1)+1, cs(2:end)];
        for j = 1:size(tinds,1)
            traces{bli}{j} = proj2one(tinds(j,1):tinds(j,2),1:10);
        end
    end
    clear Actarray

    %% Identify indices for separated trial epochs
    alcols = { t1_t2{1} ,  [6   0]   ;...
              [6     0] ,  [7   0]   ;...
    %           [7 -1000] ,  [7   0]   ;... % remove
              [7     0] ,  [8   0]   ;...
              [8     0] ,  [9   0]   ;... 
    %           [9     0] ,  [20  0]   ;... % remove
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

    %% Line up traces for each epoch and store in EPOCH
    [EPOCH] = deal(cell(length(alldays),1));
    for tb = 1:length(Project_blocks)
        bli = Project_blocks(tb); % bli = block
        for i = 1:size(alldays(bli).tt,1) % i = trial
            clc; fprintf('block: %d (%d/%d)\n',tb,i,size(alldays(bli).tt,1));
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
    % DIRREG = cellfun(@(x) [x ; nan(100*(size(x,1)<size(x,2)),10)],DIRREG,'UniformOutput',0);
    auto_EP{1} = find(cellfun(@(x) x(3),ROIs(~bad_bases))== find(sum(cellfun(@(x) x(1),alcols) == repmat([7 8],size(alcols,1),1),2)==2));
    auto_EP{2} = find(cellfun(@(x) x(3),ROIs(~bad_bases))== find(sum(cellfun(@(x) x(1),alcols) == repmat([20 20],size(alcols,1),1),2)==2));
    
    varreg = zeros(10,2);
    for di = 1:10
        varreg(di,:) = polyfit(reshape(cellfun(@(x) nanmean(x(:,di)),DIRREG),[],1),...
                                reshape(cellfun(@(x) nanstd(x(:,di)),DIRREG),[],1),1);
    end
    varfunc = @(x) x.*repmat(varreg(:,1)',size(x,1),1) + repmat(varreg(:,2)',size(x,1),1); 
    
    %% Calculate Mahalinobis distances (ALL)
    TTS = cell(length(Project_blocks),1);
    for blck  = 1:length(Project_blocks)
        TTS{blck} = alldays(blck).tt;
        M.tt{blck} = TTS{blck};
    end

    % timing stuff
    sec2min = @(tm) [floor(tm./60) floor(mod(tm,60)/10) floor(10*(mod(tm,60)/10-floor(mod(tm,60)/10)))];
    te = nan(1,sum(cellfun(@(x) size(x,1),TTS))); te(1) = 0;
    cnt = 0; 
   
    for blck = 1:length(Project_blocks)
        bli = Project_blocks(blck); % current experiment block
        T = sum(cellfun(@(x) size(x,2),EPOCH{bli}(1,:))); % time points
        P = size(alldays(bli).tt,1); % trials
        
        % initialize
        MD_t = zeros(P,T);
        MD_d = zeros(length(targs),T,P);
        MD_ep = cell(length(auto_EP),1);
%         MahrN = zeros(8,12,126,200);
        for i = 1:P
            tic; 
            cnt = cnt + 1;
            ms = sec2min((sum(cellfun(@(x) size(x,1),TTS))-cnt)*nanmean(te(1:i))); 
            clc; fprintf('%d/%d (%d:%d%d)\n',cnt,sum(cellfun(@(x) size(x,1),TTS)),ms(1),ms(2),ms(3));
%             clc; fprintf('%d/%d\n',cnt,sum(cellfun(@(x) size(x,1),TTS)));
            NTRtr = cell2mat(cellfun(@(x) x(i,:),EPOCH{bli},'UniformOutput',0))';
                
            NTR_var = varfunc(NTRtr);
            MT_boot = zeros(100,T);
            MD_boot = zeros(8,T,100);
            %%%%%
            for t = 1:T
                NTR_rnd = repmat(NTRtr(t,:),100,1) + repmat(NTR_var(t,:),100,1).*randn(100,NDIM);
                Mahr_boot = permute(reshape(cell2mat(cellfun(@(x) -sqrt(mahal(NTR_rnd,x)'),DIRREG,'UniformOutput',0)),8,100,size(roi,1)),[1 3 2]); 
                    
                [~,loc] = max(max(Mahr_boot,[],1),[],2);
                loc = squeeze(loc);
                avd = nanmean(Mahr_boot(:,mode(loc),:),3);
                
                if isnan(sum(avd))
                    MD_t(i,t) = NaN;
                    MD_d(:,t,i) = nan(8,1);
                else
                    MD_t(i,t) = mode(loc);
                    MD_d(:,t,i) = avd;
                end
            end
                
%                 
%             for tboot = 1:100
%                 NTR_rnd = NTRtr + NTR_var.*randn(size(NTRtr));
%             
%                 Mahr_boot = permute(reshape(cell2mat(cellfun(@(x) -sqrt(mahal(NTR_rnd,x)'),DIRREG,'UniformOutput',0)),8,T,[]),[1 3 2]);
% %                 Mahr = permute(reshape(cell2mat(cellfun(@(x) -sqrt(mahal(NTR{i},x)'),DIRREG,'UniformOutput',0)),8,T,[]),[1 3 2]);
%                 for t = 1:T
%                     Maldist_boot = Mahr_boot(:,:,t);
%                     maxoverdir = max(Maldist_boot,[],1);
%                     diratmax = nanmean(Maldist_boot(:,maxoverdir==max(maxoverdir)),2);
%                     
%                     if isempty(diratmax)
%                         MT_boot(tboot,t) = nan;
%                         MD_boot(:,t,tboot) = nan(8,1);
%                     else
%                         [~,maxtm] = max(maxoverdir);
%                         MT_boot(tboot,t) = roicents(maxtm);
%                         MD_boot(:,t,tboot) = circshift(diratmax,5-targids{bli}(i));
%                     end
%                 end
%             end
%             MD_t(i,:) = nanmean(MT_boot);
%             MD_d(:,:,i) = nanmean(MD_boot,3);
%                
%                 
%                 for q = 1:length(auto_EP)
%                     if isnan(sum(Maldists(:,auto_EP{q})))
%                         MD_ep{q}(:,t,i) = nan(8,1);
%                     else
%                         MD_ep{q}(:,t,i) = circshift(Maldists(:,auto_EP{q}),5-targids{bli}(i));
%                     end
%                 end

            

            te(cnt) = toc;
        end
        M.([lower(area) '_noise2']){blck}.d = MD_d;
        M.([lower(area) '_noise2']){blck}.t = MD_t;
        
    end

end
M.meta.('filter_size') = filtsizes;
M.meta.('time_bin_size') = dsamplesize;
M.meta.('filter_type') = filttype;
clc; fprintf('\n----------\n-- Done --\n----------\n');