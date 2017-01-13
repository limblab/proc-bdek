%% Parameters
areas = {'PMd','M1'};%,'M1'};%'M1'};%,'M1'};%,'M1'};%,'M1'};%,'M1'};%,'M1'};%{'PMd','M1'};
NBOOT = 1;
do_square_root_transform = false;
include_aborted = false; % 'aborted';

NDIM = 10;
Baseline_block = 1;
dsamplesize = 25; % In milliseconds
filtsize = 150; % In milliseconds
filttype = 'causal';
align_col = 19; % 17 for chosen, 13 for correct

T1_T2 = { [6  -200] , [6   0] ;...
          [6     0] , [7   0] ;...
          [7     0] , [8   0] ;...
          [8     0] , [9   0] ;...
          [9     0] , [3   0] };    
      
% T1_T2 = { [5  -200] , [5   0] ;...
%           [5     0] , [6   0] ;...
%           [6     0] , [7   0] }; 
%       
useforcomp = 2:5;

%% Preprocess trial table
if ~exist('TT','var')
    TT = cell(1,20);
    for i = 1:length(alldays)
        TT{i} = alldays(i).tt;
        if include_aborted && ~isempty(alldays(i).ttA)
            TT{i+length(alldays)} = alldays(i).ttA;
        end
    end
    TT(cellfun(@isempty,TT)) = [];
    for i = 1:length(TT)
        % Remove trials missing timestamps or target direction
    %     TT{i}(isnan(sum(TT{i}(:,unique(reshape(cell2mat(cellfun(@(x) x(1),T1_T2,'Uni',0)),[],1))),2)),:) = [];
        TT{i}(isnan(TT{i}(:,align_col)),:) = [];
    end
end
targids = cellfun(@(x) mod(round(mod(x(:,align_col)+2*pi,2*pi)./(pi/4)),8)+1,TT,'Uni',0);
dirinds = cell(length(unique(targids{Baseline_block})),1);
for i = unique(targids{Baseline_block})'
    dirinds{i} = find(targids{Baseline_block}==i);
end
trialkey = zeros(length(targids{Baseline_block}),2);
for i= 1:length(targids{Baseline_block})
    trialkey(i,:) = [targids{Baseline_block}(i), find(dirinds{targids{Baseline_block}(i)} == i)];
end

%% Set up activity cell
Actarray = [];
for ba_loop = 1:length(areas) % Brain area
    Actarray.(areas{ba_loop}) = cell(length(TT),1);
    for trial_block = 1:length(TT) % trial block
        Actarray.(areas{ba_loop}){trial_block} = cell(size(TT{trial_block},1),size(T1_T2,1));
        for epoch = 1:size(T1_T2,1) % trial epoch (planning, execution, etc.)  
            
            clc; fprintf('area: %s\nblock: %d/%d\nepoch: %d/%d\n',areas{ba_loop},trial_block,length(TT),epoch,size(T1_T2,1));
            Actarray.(areas{ba_loop}){trial_block}(:,epoch) = ...
                    trial_raster_withsmoothing(eval(sprintf('alldays(1).%s_units',areas{ba_loop})),... % brain area
                                                TT{trial_block},T1_T2{epoch,1},T1_T2{epoch,2},...       % Trial table info
                                                filtsize,strcmp(filttype,'causal'),dsamplesize);        % Filter params        
        end
    end
end
Actarray.meta.('filter_size') = filtsize;
Actarray.meta.('time_bin_size') = dsamplesize;
Actarray.meta.('filter_type') = filttype;

%% do PCA and project activity
[PCax,eigenvals,mu] = deal(cell(length(areas),1));
Actproj = [];
for ba_loop = 1:length(areas)
    [PCax{ba_loop},~,eigenvals{ba_loop},~,~,mu{ba_loop}] = pca(cell2mat(reshape(vertcat(Actarray.(areas{ba_loop}){:}),1,[]))'); 
    for epoch = 1:length(Actarray.(areas{ba_loop})) %% trial block
        Actproj.(areas{ba_loop}){epoch} = cellfun(@(x) ((x'-repmat(mu{ba_loop},size(x,2),1))*PCax{ba_loop}(:,1:NDIM))',Actarray.(areas{ba_loop}){epoch},'Uni',0);
    end
end

%% Set up directional sets
DR = []; DRbl = [];
for ba_loop = 1:length(areas)
    for i = 1:length(dirinds)
        DR.(areas{ba_loop}){i,1} = cell2mat(reshape(Actproj.(areas{ba_loop}){Baseline_block}(targids{Baseline_block}==i,useforcomp),1,[]))';
        DR.(areas{ba_loop}){i}(isnan(sum(DR.(areas{ba_loop}){i},2)),:) = [];
        for j = 1:length(dirinds{i})
            DRbl.(areas{ba_loop}){i,1}{j,:} = cell2mat(Actproj.(areas{ba_loop}){Baseline_block}(dirinds{i}(j),useforcomp));
            DRbl.(areas{ba_loop}){i}{j}(:,isnan(sum(DRbl.(areas{ba_loop}){i}{j},1))) = [];
        end
    end
end

%% Calculate distances
D2base = [];
for ba_loop = 1:length(areas)

    for trial_block = 1:length(Actproj.(areas{ba_loop}))

        if trial_block == Baseline_block
            for j = 1:size(Actproj.(areas{ba_loop}){Baseline_block},1)
                clc; fprintf('area: %s\nblock: %d/%d\nleave-one-out for baseline block: %d/%d\n',...
                             areas{ba_loop},trial_block,length(Actproj.(areas{ba_loop})),j,size(Actproj.(areas{ba_loop}){Baseline_block},1))
                tempDR = DRbl.(areas{ba_loop});
                tempDR{trialkey(j,1)}{trialkey(j,2)} = [];
                
                for i = unique(targids{Baseline_block}')
                    D2base.(areas{ba_loop}){trial_block}{i,:}(j,:) = cellfun(@(x) mahal(x',cell2mat(tempDR{i}')')',Actproj.(areas{ba_loop}){trial_block}(j,:),'Uni',0);
                end
            end
 
        else
            clc; fprintf('area: %s\nblock: %d/%d\n',...
                             areas{ba_loop},trial_block,length(Actproj.(areas{ba_loop})));
            for i = unique(targids{Baseline_block}')
                D2base.(areas{ba_loop}){trial_block}{i,:} = cellfun(@(x) mahal(x',DR.(areas{ba_loop}){i})',Actproj.(areas{ba_loop}){trial_block},'Uni',0);
            end
        end
    end
end
clear tempDR     

%% Trial-by-trial distances
Dtrial = [];
for ba_loop = 1:length(areas)

    for trial_block = 1:length(Actproj.(areas{ba_loop}))
        
        for trial = 1:size(Actproj.(areas{ba_loop}){trial_block},1)
            clc; fprintf('%s\nblock: %d/%d\ntrial: %d/%d\n',areas{ba_loop},trial_block,...
                length(Actproj.(areas{ba_loop})),trial,size(Actproj.(areas{ba_loop}){trial_block},1));
            
            a = cellfun(@(x) x(trial,:),D2base.(areas{ba_loop}){trial_block},'Uni',0);
            asizes = cellfun(@(x) size(x,2),a{1});
            b = 1./sqrt(circshift(cell2mat(vertcat(a{:})),5-targids{trial_block}(trial),1));
            
            Dtrial.(areas{ba_loop}){trial_block}(trial,:) =  mat2cell(b,size(b,1),asizes);
        end

    end
end

%% Padded distances
Dpad = []; Dcomb = [];
epsizes = max(cellfun(@(x) size(x,2), vertcat(Dtrial.(areas{1}){:})));
LENS = epsizes;
CLENS = [1 cumsum(epsizes)];
for ba_loop = 1:length(areas)
    for trial_block = 1:length(Actproj.(areas{ba_loop}))
        for trial = 1:size(Actproj.(areas{ba_loop}){trial_block},1)
            clc; fprintf('%s\nblock: %d/%d\ntrial: %d/%d\n',areas{ba_loop},trial_block,...
                length(Actproj.(areas{ba_loop})),trial,size(Actproj.(areas{ba_loop}){trial_block},1));
            for ep = 1:length(Dtrial.(areas{ba_loop}){trial_block}(trial,:))
                Dpad.(areas{ba_loop}){trial_block}{1,ep}(:,:,trial) = [Dtrial.(areas{ba_loop}){trial_block}{trial,ep}, ...
                                                               NaN(8,epsizes(ep)-size(Dtrial.(areas{ba_loop}){trial_block}{trial,ep},2))];
            end
        end
        Dcomb.(areas{ba_loop}){trial_block} = cat(2,Dpad.(areas{ba_loop}){trial_block}{:});

    end
end
Dcomb.bounds = [1 cumsum(epsizes)];

%% Padded to median
Dmed = []; DmedT = [];
epsizesm = median(cellfun(@(x) size(x,2), vertcat(Dtrial.(areas{1}){:})));
LENSmed = epsizesm;
CLENSmed = [1 cumsum(epsizesm)];
for ba_loop = 1:length(areas)
    for trial_block = 1:length(Actproj.(areas{ba_loop})) 
        for trial = 1:size(Actproj.(areas{ba_loop}){trial_block},1)
            clc; fprintf('%s\nblock: %d/%d\ntrial: %d/%d\n',areas{ba_loop},trial_block,...
                length(Actproj.(areas{ba_loop})),trial,size(Actproj.(areas{ba_loop}){trial_block},1));
            for ep = 1:length(Dtrial.(areas{ba_loop}){trial_block}(trial,:))
                if size(Dtrial.(areas{ba_loop}){trial_block}{trial,ep},2) < epsizesm(ep)
                    DmedT.(areas{ba_loop}){trial_block}{1,ep}(:,:,trial) = [Dtrial.(areas{ba_loop}){trial_block}{trial,ep}, ...
                                                                   NaN(8,epsizesm(ep)-size(Dtrial.(areas{ba_loop}){trial_block}{trial,ep},2))];
                else
                    DmedT.(areas{ba_loop}){trial_block}{1,ep}(:,:,trial) = Dtrial.(areas{ba_loop}){trial_block}{trial,ep}(:,1:epsizesm(ep));
                end
            end
        end
        Dmed.(areas{ba_loop}){trial_block} = cat(2,DmedT.(areas{ba_loop}){trial_block}{:});
    end
end
Dmed.bounds = [1 cumsum(epsizesm)];
clear DmedT

%% Compile d metrics into single variable 'D'
D = [];
D.full = Dtrial;
D.med = Dmed;
D.tt = TT;
D.targids = targids;
D.meta = Actarray.meta;

% M.CLENS = CLENS;
% M.LENS = LENS;

clc; fprintf('\n----------\n-- Done --\n----------\n');

%% AB for each region
AB = [];
ABfull = [];
for ba_loop = 1:length(areas)
    for trial_block = 1:length(Dtrial.(areas{ba_loop})) 
        ABS = cellfun(@(x) nanmean(x([5 1],:)./repmat(nanmean(x([3 7],:),1),2,1),2)',Dtrial.(areas{ba_loop}){trial_block},'Uni',0);
        AB.(areas{ba_loop}){trial_block} = mat2cell(cell2mat(ABS),size(ABS,1),2*ones(size(ABS,2),1));
        
        ABSf = cellfun(@(x) (x([5 1],:)./repmat(nanmean(x([3 7],:),1),2,1))',Dtrial.(areas{ba_loop}){trial_block},'Uni',0);
        for i = 1:size(ABSf,2)
            ABfull.(areas{ba_loop}){trial_block}{1,i} = cell2mat(ABSf(:,i));
        end
    end
end

            
            

            
            