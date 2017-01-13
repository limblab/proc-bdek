function[rasters, fig_hand, graph_objs, kin] = PSTH_heat(bdf,train,trial_table,t1_t2,WindowBin,DirBin,alignment)
%%%

if strcmp(alignment,'target')
    aligntype = 6; % 6 for target, 7 for go cue
elseif strcmp(alignment,'go')
    aligntype = 7; 
elseif strcmp(alignment,'reward');
    aligntype = 8;
end

t1 = t1_t2(1); t2 = t1_t2(2);

% Pick out rows from trial table for low and high uncertainty
lowinds = find(trial_table(:,3)==max(trial_table(:,3)));
highinds = find(trial_table(:,3)==min(trial_table(:,3)));
lowvartrials = [lowinds trial_table(lowinds,:)];
highvartrials = [highinds trial_table(highinds,:)];

winedges = 1:WindowBin:round(1000*(t2-t1))+1;
wincents = t1*1000+WindowBin/2:WindowBin:t2*1000-WindowBin/2;

rasters{1} = zeros(length(lowvartrials),round(1000*(t2-t1)));
rasters{2} = zeros(length(highvartrials),round(1000*(t2-t1)));

kin{1}.posx  = zeros(length(lowvartrials),round(1000*(t2-t1)));
kin{1}.posy  = zeros(length(lowvartrials),round(1000*(t2-t1)));
kin{1}.velx  = zeros(length(lowvartrials),round(1000*(t2-t1)));
kin{1}.vely  = zeros(length(lowvartrials),round(1000*(t2-t1)));
kin{1}.speed = zeros(length(lowvartrials),round(1000*(t2-t1)));

kin{2}.posx  = zeros(length(highvartrials),round(1000*(t2-t1)));
kin{2}.posy  = zeros(length(highvartrials),round(1000*(t2-t1)));
kin{2}.velx  = zeros(length(highvartrials),round(1000*(t2-t1)));
kin{2}.vely  = zeros(length(highvartrials),round(1000*(t2-t1)));
kin{2}.speed = zeros(length(highvartrials),round(1000*(t2-t1)));

%% Fill out Rasters
lcount = 1;
hcount = 1;
for i = 1:length(trial_table)
    
    % Get Time
    timestart = trial_table(i,aligntype-1)+t1;
    feedback = trial_table(i,3)==min(trial_table(:,3));
   
    bdfstart = find(bdf.pos(:,1)<timestart,1,'last');
    bdfend = bdfstart + round(1000*(t2-t1))-1;
    
    %%% Kinematics
    pos = bdf.pos(bdfstart:bdfend,2:3)';
    vel = bdf.vel(bdfstart:bdfend,2:3)';
    speed = sqrt(vel(1,:).^2 + vel(2,:).^2);
    
    % Align spikes to start and get rid of those not in region
    aligned_ts = round(1000*(train - timestart));
    aligned_ts(aligned_ts<=0 | aligned_ts>=((t2-t1)*1000)) = [];
    
    if feedback == 0
        rasters{1}(lcount,aligned_ts) = 1;
        
        kin{1}.posx(lcount,:) = pos(1,:);
        kin{1}.posy(lcount,:) = pos(2,:);
        kin{1}.velx(lcount,:) = vel(1,:);
        kin{1}.vely(lcount,:) = vel(2,:);
        kin{1}.speed(lcount,:)= speed;
        
        lcount = lcount + 1;
        
    elseif feedback == 1
        rasters{2}(hcount,aligned_ts) = 1;
        
        kin{2}.posx(hcount,:) = pos(1,:);
        kin{2}.posy(hcount,:) = pos(2,:);
        kin{2}.velx(hcount,:) = vel(1,:);
        kin{2}.vely(hcount,:) = vel(2,:);
        kin{2}.speed(hcount,:)= speed;
        
        hcount = hcount + 1;
    
    else
        warning('Check Trial Table');
    end
    
end

%%
dispersions = unique(trial_table(:,3));
disps = cell(length(dispersions),1);
for i = 1:length(dispersions)
    disps{i} = sprintf('k = %d',dispersions(i));
end
%% Fill out bins of PSTH
Lwin = zeros(length(lowvartrials),length(winedges)-1);
Hwin = zeros(length(highvartrials),length(winedges)-1);
for i = 1:length(winedges)-1
    Lwin(:,i) = sum(rasters{1}(:,winedges(i):winedges(i+1)-1),2);
    Hwin(:,i) = sum(rasters{2}(:,winedges(i):winedges(i+1)-1),2);
end
Lwin = 1000*Lwin./WindowBin;
Hwin = 1000*Hwin./WindowBin;

L_dirs = trial_table(lowinds,10);
H_dirs = trial_table(highinds,10);

bins = 0:DirBin*pi/180:2*pi; bins(end) = 10;

binned_L = zeros(length(bins)-1,size(Lwin,2));
binned_H = zeros(length(bins)-1,size(Hwin,2));
for i = 1:length(bins)-1
    
    binned_L(i,:) = nanmean(Lwin(L_dirs>bins(i) & L_dirs<bins(i+1),:),1);
    binned_H(i,:) = nanmean(Hwin(H_dirs>bins(i) & H_dirs<bins(i+1),:),1);
    
end


bincounts_L = mean(Lwin,1);
bincount_SE_L = std(Lwin,[],1)./sqrt(length(lowvartrials));

bincounts_H = mean(Hwin,1);
bincount_SE_H = std(Hwin,[],1)./sqrt(length(highvartrials));

%%
fig_hand = figure; hold on;

if length(unique(trial_table(:,3)))>1
    plot(wincents, bincounts_H , 'r.-');    
    plot(wincents, bincounts_L , 'b.-');
    
    legend(disps);
    
    patch([wincents fliplr(wincents)],...
        [bincounts_H - 1.96*bincount_SE_H, fliplr(bincounts_H + 1.96*bincount_SE_H)],'r','FaceAlpha',0.4,'EdgeAlpha',0.1);
    patch([wincents fliplr(wincents)],...
        [bincounts_L - 1.96*bincount_SE_L, fliplr(bincounts_L + 1.96*bincount_SE_L)],'b','FaceAlpha',0.4,'EdgeAlpha',0.1);
else
    plot(wincents, bincounts_H , 'k.-');
    patch([wincents fliplr(wincents)],...
        [bincounts_H - 1.96*bincount_SE_H, fliplr(bincounts_H + 1.96*bincount_SE_H)],'k','FaceAlpha',0.4,'EdgeAlpha',0.1);
end

max_ylim = ylim;
plot([0 0],[0 200],'k--');

graph_objs.legend = disps;
graph_objs.axis = [min(wincents) max(wincents) 0 max_ylim(2)+2];
graph_objs.xlabel = sprintf('Time from %s (ms)',alignment);
graph_objs.ylabel = 'Firing Rate (/sec)';
graph_objs.prior_mean = mean_angle(trial_table(:,2),'rads');
graph_objs.prior_kappa = calc_disp(trial_table(:,2));

if graph_objs.prior_kappa < 0.25
    graph_objs.prior_kappa = 0;
    graph_objs.prior_mean = 0;
end

axis(graph_objs.axis);
xlabel(graph_objs.xlabel,'FontSize',14);
ylabel(graph_objs.ylabel,'FontSize',14);

if sum(sum(kin{1}.speed)) ==0 
    kin{1} = kin{2};
end


