function[OUTPUTS,Lout,Hout] = raster_UNT_circ(bdf,train,trial_table,t1,t2,alignment,do_plot,position_rank)
%%%
% OUTPUTS = [meanFR, final_pos, target_pos, maxspeed, maxspeed_time,
% direction, UNC];

if strcmp(alignment,'target')
    aligntype = 6; % 6 for target, 7 for go cue
elseif strcmp(alignment,'go')
    aligntype = 7;
else 
    fprintf('choose ''target'' or ''go''\n');
end


if do_plot==1; Plotraster = 1; else Plotraster = 0; end% 1 to plot raster

Plotcont = 0; %1 to plot continuous firing rate estimates
subcont = 0;

Plotspeeds = 0; % 1 to plot average speed profiles on top

WindowBin = 0; 
Endpoint_rank = 1;
PosBinned = 0;
BlobPlot = 0;

if strcmp(position_rank,'endpoint')
    rank_anchor = 1; 
elseif strcmp(position_rank,'target')
    rank_anchor = 2;
elseif strcmp(position_rank,'none')
    Endpoint_rank = 0;
    rank_anchor = 0;
end


                 % 1: end position
                 % 2: target position
                 % 3: centroid position NOT AVAILABLE YET
                 % 4: angle at max speed 
                 % 0: ONLY USE TO GET FULL OUTPUT MATRIX
emprow = 0;
rastershift = 0;    
contParam = 100;
xs = [-10 2 12]; % Range for endpoint/FR plots
%xs = [0 1.6149 pi];
Errorbarbinplot = 0;

% Pick out rows from trial table for low and high uncertainty
lowinds = find(trial_table(:,3)==max(trial_table(:,3)));
highinds = find(trial_table(:,3)==min(trial_table(:,3)));
lowvartrials = [lowinds trial_table(lowinds,:)];
highvartrials = [highinds trial_table(highinds,:)];

LL = length(lowvartrials);
LH = length(highvartrials);

if WindowBin > 0
    winedges = floor(1:WindowBin:round(1000*(t2-t1))+2000);
    wincents = 0.5*WindowBin:WindowBin:winedges(end)-0.5*WindowBin-1;
end

%% Low Var Trials
if Plotraster == 1
	figure; hold on;
end
Lraster = nan(LL,round(1000*(t2-t1))); % Raster with NaN
Lraster_ext = nan(LL,round(1000*(t2-t1))+2000);
Lbin = zeros(LL,round(1000*(t2-t1))); % Raster with Zeros
Lcont = zeros(LL,round(1000*(t2-t1))); % Gauss Cont raster
Lspeed = zeros(size(Lcont)); % Speeds
Lpos_spikes = nan(LL,2);
L_ISI = cell(LL,1);
L_LP = nan(LL,3);

for i = 1:LL
    
    % Get Time of start
    timestart = lowvartrials(i,aligntype)+t1;
    bdfstart = find(bdf.pos(:,1)<timestart,1,'last');
    bdfgo = find(bdf.pos(:,1)<lowvartrials(i,7),1,'last');
    bdfend = bdfstart + 1000*(t2-t1)-1;
    bdfrew = find(bdf.pos(:,1)>lowvartrials(i,8),1,'first');
    
    % Get Endpoint/target pos/initial angle
    if rank_anchor == 1
        final_pos = atan2(bdf.pos(bdfrew,3),bdf.pos(bdfrew,2));
    elseif rank_anchor == 2
        final_pos = mod(lowvartrials(i,3),2*pi);
    elseif rank_anchor == 3
        final_pos = 0;
    elseif rank_anchor == 4
        speed = sqrt((bdf.vel(bdfgo:bdfrew,2)').^2 + (bdf.vel(bdfgo:bdfrew,3)').^2);
        maxspeedind = find(speed==max(speed));
        maxvelvec = bdf.vel(bdfgo+maxspeedind-1,[2 3]);
        direction = atan2(maxvelvec(2),maxvelvec(1));
        final_pos = direction;
    elseif rank_anchor == 0
        L_final_pos(i) = bdf.pos(bdfrew,2);
        L_target_pos(i) = lowvartrials(i,3);
        
        speed = sqrt((bdf.vel(bdfgo:bdfrew,2)').^2 + (bdf.vel(bdfgo:bdfrew,3)').^2);
        L_maxspeed(i) = max(speed);
        maxspeedind = find(speed==max(speed));
        L_maxspeed_time(i) = bdfgo-bdfstart+maxspeedind-1;
        maxvelvec = bdf.vel(bdfgo+maxspeedind-1,[2 3]);
        L_direction(i) = atan2(maxvelvec(2),maxvelvec(1));
        
        final_pos = L_direction(i);
    end
    
    % Find Speed Profile of trial
    Lspeed(i,:) = sqrt((bdf.vel(bdfstart:bdfend,2)').^2 + (bdf.vel(bdfstart:bdfend,3)').^2);
    
    % Align spikes to start and get rid of those not in region
    aligned_ts = round(1000*(train - timestart));
    aligned_ts(aligned_ts<=0 | aligned_ts>=((t2-t1)*1000)) = [];
    
    aligned_ts_exp = round(1000*(train -(timestart-1)));
    aligned_ts_exp(aligned_ts_exp<=0 | aligned_ts_exp>=((t2-t1+2)*1000)) = [];
    
    % Fill out rasters
    Lraster(i,aligned_ts) = 1;
    Lbin(i,aligned_ts) = 1;
    
    tempraster = zeros(1,1000*(t2-t1+2));
    tempraster(aligned_ts_exp) = 1;
    Lraster_ext(i,:) = tempraster;
    
    % Convert binned spikes to continuous (Gaussian) firing rate
    tempcont = train2cont(tempraster,contParam);
    tempcont(1:1000) = []; tempcont(end-999:end) = [];
    Lcont(i,:) = tempcont;
        
    Lpos_spikes(i,:) = [final_pos sum(Lbin(i,:))];
    L_ISI{i} = diff(aligned_ts);
    
    L_LP(i,:) = [final_pos find(Lcont(i,:)==max(Lcont(i,:)),1,'first') max(Lcont(i,:))];
    
    if Plotraster == 1
        if Endpoint_rank == 0
            plot((1000*t1):1000*t2-1,i*Lraster(i,:),'b.');
            if Lpos_spikes(i,2)==0 && emprow == 1
                plot((1000*t1):1000*t2-1,i*ones(1,length(t2-t1)),'bx','MarkerSize',1);
            end    
            
        else
            plot((1000*t1):1000*t2-1,final_pos*Lraster(i,:),'b.');
            if Lpos_spikes(i,2)==0 && emprow == 1
                plot((1000*t1):1000*t2-1,final_pos*ones(1,length(t2-t1)),'bx','MarkerSize',1);
            end   
        end
    end
    
    
end

Lavspeed = mean(Lspeed,1); Lstdspeed = std(Lspeed,1); maxLavs = max(Lavspeed);
if Plotspeeds==1 && Plotraster==1
    plot((1000*t1):1000*t2-1,(LL/maxLavs).*Lavspeed,'k');
end

%% High Var Trials
Hraster = nan(length(highvartrials),round(1000*(t2-t1))); % Raster with NaN
Hraster_ext = nan(LH,round(1000*(t2-t1))+2000);
Hbin = zeros(LH,round(1000*(t2-t1))); % Raster with Zeros
Hcont = zeros(LH,round(1000*(t2-t1))); % Gauss Cont raster
Hspeed = zeros(size(Hcont)); % Speeds
Hpos_spikes = nan(LH,2);
H_ISI = cell(LH,1);
H_LP = nan(LH,3);
for i = 1:LH
    
    % Get Time of start
    timestart = highvartrials(i,aligntype)+t1;
    bdfstart = find(bdf.pos(:,1)<timestart,1,'last');
    bdfgo = find(bdf.pos(:,1)<highvartrials(i,7),1,'last');
    bdfend = bdfstart+round(1000*(t2-t1)-1);
    bdfrew = find(bdf.pos(:,1)>highvartrials(i,8),1,'first');
    
    % Get Endpoint/target pos/initial angle
    if rank_anchor == 1
        final_pos = atan2(bdf.pos(bdfrew,3),bdf.pos(bdfrew,2));
    elseif rank_anchor == 2
        final_pos = mod(highvartrials(i,3),2*pi);
    elseif rank_anchor == 3
        final_pos = 0;
    elseif rank_anchor == 4
        speed = sqrt((bdf.vel(bdfgo:bdfrew,2)').^2 + (bdf.vel(bdfgo:bdfrew,3)').^2);
        maxspeedind = find(speed==max(speed));
        maxvelvec = bdf.vel(bdfgo+maxspeedind-1,[2 3]);
        direction = atan2(maxvelvec(2),maxvelvec(1));
        final_pos = direction;
   elseif rank_anchor == 0
        H_final_pos(i) = bdf.pos(bdfrew,2);
        H_target_pos(i) = highvartrials(i,3);
        
        speed = sqrt((bdf.vel(bdfgo:bdfrew,2)').^2 + (bdf.vel(bdfgo:bdfrew,3)').^2);
        H_maxspeed(i) = max(speed);
        maxspeedind = find(speed==max(speed));
        H_maxspeed_time(i) = bdfgo-bdfstart+maxspeedind-1;
        maxvelvec = bdf.vel(bdfgo+maxspeedind-1,[2 3]);
        H_direction(i) = atan2(maxvelvec(2),maxvelvec(1));
        final_pos = H_direction(i);
    end
    
    % Find Speed Profile of trial
    Hspeed(i,:) = sqrt((bdf.vel(bdfstart:bdfend,2)').^2 + (bdf.vel(bdfstart:bdfend,3)').^2);
    
    % Align spikes to start and get rid of those not in region
    aligned_ts = round(1000*(train - timestart));
    aligned_ts(aligned_ts<=0 | aligned_ts>=((t2-t1)*1000)) = [];
    aligned_ts_exp = round(1000*(train - (timestart - 1)));
    aligned_ts_exp(aligned_ts_exp<=0 | aligned_ts_exp>=((t2-t1+2)*1000)) = [];
    
    % Fill out rasters
    Hraster(i,aligned_ts) = 1;
    Hbin(i,aligned_ts) = 1;
    
    tempraster = zeros(1,1000*(t2-t1+2));
    tempraster(aligned_ts_exp) = 1;
    Hraster_ext(i,:) = tempraster;
    
    % Convert binned spikes to continuous firing rate
    tempcont = train2cont(tempraster,contParam);
    tempcont(1:1000) = []; tempcont(end-999:end) = [];
    Hcont(i,:) = tempcont;
    
    Hpos_spikes(i,:) = [final_pos sum(Hbin(i,:))];
    H_ISI{i} = diff(aligned_ts);
    
    H_LP(i,:) = [final_pos find(Hcont(i,:)==max(Hcont(i,:)),1,'first') max(Hcont(i,:))];

    
    if Plotraster == 1 && (sum(lowinds(1)==highinds(1))==0)
        if Endpoint_rank == 0
            plot((1000*t1):1000*t2-1,LL+i*Hraster(i,:),'r.');
            if Hpos_spikes(i,2)==0 && emprow == 1
                plot((1000*t1):1000*t2-1,LL+i*ones(1,length(t2-t1)),'rx','MarkerSize',1);
            end  
        else
            plot((1000*t1):1000*t2-1,rastershift+final_pos*Hraster(i,:),'r.');
            if Hpos_spikes(i,2)==0 && emprow == 1
                plot((1000*t1):1000*t2-1,rastershift+final_pos*ones(1,length(t2-t1)),'rx','MarkerSize',1);
            end   
        end
    end
end


Havspeed = mean(Hspeed,1); Hstdspeed = std(Hspeed,1); maxHavs = max(Havspeed);
if Plotspeeds ==1 && Plotraster==1
    plot((1000*t1):1000*t2-1,LL+(LH/maxHavs).*Havspeed,'k');
end

%% Setup output if rank anchor == 0

% L_meancont = mean(Lcont,2);
% H_meancont = mean(Hcont,2);

L_meancont = Lpos_spikes(:,2);
H_meancont = Hpos_spikes(:,2);

if Plotraster == 1
    Lout = Lraster;
    Hout = Hraster;
elseif Plotcont == 1
    Lout = Lcont;
    Hout = Hcont;
else
    Lout = Lraster;
    Hout = Hraster;
end

if rank_anchor == 0

    OUTPUTS.x = zeros(LL+LH,6);
    OUTPUTS.y(lowinds,:) = L_meancont;
    OUTPUTS.y(highinds,:) = H_meancont;
    
    OUTPUTS.x(lowinds,:) = [L_final_pos', L_target_pos', ...
                          L_maxspeed', L_maxspeed_time', L_direction', ones(LL,1)];
                      
    OUTPUTS.x(highinds,:) = [H_final_pos', H_target_pos', ...
                          H_maxspeed', H_maxspeed_time', H_direction', 2*ones(LH,1)];
                      
elseif rank_anchor ~= 0
    
    OUTPUTS.y = zeros(LL+LH,1);
    
    OUTPUTS.y(lowinds) = L_meancont;
    OUTPUTS.y(highinds) = H_meancont;
    
end

%% Plot Continuous FR and/or Speeds

alignopts = cell(2,1);
alignopts{1} = 'Target On';
alignopts{2} = 'Go Cue';

if Plotcont == 1 && Plotspeeds==1 && subcont==1

    figure; subplot(2,1,1); hold on;
    
    % FR mean and std
    Lcont_mean = mean(Lcont,1);
    Lcont_std = std(Lcont,1);
    Hcont_mean = mean(Hcont,1);
    Hcont_std = std(Hcont,1);
  
    % Do Plots of FR
    plot((1000*t1):1000*t2-1,Lcont_mean,'b','LineWidth',3);
    plot((1000*t1):1000*t2-1,Lcont_mean+Lcont_std,'b','LineWidth',1);
    %plot((1000*t1):1000*t2-1,Lcont_mean-Lcont_std,'b','LineWidth',1);
    
    plot((1000*t1):1000*t2-1,Hcont_mean,'r','LineWidth',3);
    plot((1000*t1):1000*t2-1,Hcont_mean+Hcont_std,'r','LineWidth',1);
    %plot((1000*t1):1000*t2-1,Hcont_mean-Hcont_std,'r','LineWidth',1);
    
    title('Neural Activity','FontSize',16);
    xlabel(sprintf('Time from %s (ms)',alignopts{aligntype-5}),'FontSize',14);
    ylabel('Firing Rate (spikes/sec)','FontSize',14);
    
    % Do Plots of Speed
    subplot(2,1,2); hold on;
    plot((1000*t1):1000*t2-1,Lavspeed,'b','LineWidth',3);
    plot((1000*t1):1000*t2-1,Lavspeed+Lstdspeed,'b','LineWidth',1);
    plot((1000*t1):1000*t2-1,Lavspeed-Lstdspeed,'b','LineWidth',1);
    
    plot((1000*t1):1000*t2-1,Havspeed,'r','LineWidth',3);
    plot((1000*t1):1000*t2-1,Havspeed+Hstdspeed,'r','LineWidth',1);
    plot((1000*t1):1000*t2-1,Havspeed-Hstdspeed,'r','LineWidth',1);

    title('Hand Speed','FontSize',16);
    xlabel(sprintf('Time from %s (ms)',alignopts{aligntype-5}),'FontSize',14);
    ylabel('Speed (cm/s)','FontSize',14);
    
elseif Plotcont == 1 && Plotspeeds ~= 1 && subcont==1

    figure; 
    
    Lcont_mean = mean(Lcont,1);
    Lcont_std = std(Lcont,1);
    Hcont_mean = mean(Hcont,1);
    Hcont_std = std(Hcont,1);
  
    subplot(2,1,1); hold on;
    plot((1000*t1):1000*t2-1,Lcont_mean,'b','LineWidth',3);
    plot((1000*t1):1000*t2-1,Lcont_mean+Lcont_std,'b','LineWidth',1);
    %plot((1000*t1):1000*t2-1,Lcont_mean-Lcont_std,'b','LineWidth',1);
    
    plot((1000*t1):1000*t2-1,Hcont_mean,'r','LineWidth',3);
    plot((1000*t1):1000*t2-1,Hcont_mean+Hcont_std,'r','LineWidth',1);
    %plot((1000*t1):1000*t2-1,Hcont_mean-Hcont_std,'r','LineWidth',1);
    
    title('Neural Activity','FontSize',16);
    xlabel(sprintf('Time from %s (ms)',alignopts{aligntype-5}),'FontSize',14);
    ylabel('Firing Rate (spikes/sec)','FontSize',14);
    
    subplot(2,1,2); hold on;
    plot((1000*t1):1000*t2-1,Lcont_std.^2,'b','LineWidth',3);
    plot((1000*t1):1000*t2-1,Hcont_std.^2,'r','LineWidth',3);

    title('FR Variance','FontSize',16);
    xlabel(sprintf('Time from %s (ms)',alignopts{aligntype-5}),'FontSize',14);
    ylabel('Firing Rate variance','FontSize',14);
    
elseif Plotcont ~= 1 && Plotspeeds == 1 && subcont==1
   
    figure; hold on;
    
    plot((1000*t1):1000*t2-1,Lavspeed,'b','LineWidth',3);
    plot((1000*t1):1000*t2-1,Lavspeed+Lstdspeed,'b','LineWidth',1);
    plot((1000*t1):1000*t2-1,Lavspeed-Lstdspeed,'b','LineWidth',1);
    
    plot((1000*t1):1000*t2-1,Havspeed,'r','LineWidth',3);
    plot((1000*t1):1000*t2-1,Havspeed+Hstdspeed,'r','LineWidth',1);
    plot((1000*t1):1000*t2-1,Havspeed-Hstdspeed,'r','LineWidth',1);
    
    title('Hand Speed','FontSize',16);
    xlabel(sprintf('Time from %s (ms)',alignopts{aligntype-5}),'FontSize',14);
    ylabel('Speed (cm/s)','FontSize',14);
    
elseif Plotcont == 1 && subcont==0
    figure; hold on;
    
    % FR mean and std
    Lcont_mean = mean(Lcont,1);
    Lcont_std = std(Lcont,1);
    Hcont_mean = mean(Hcont,1);
    Hcont_std = std(Hcont,1);
  
    % Do Plots of FR
    plot((1000*t1):1000*t2-1,Lcont_mean,'b','LineWidth',3);
    %plot((1000*t1):1000*t2-1,Lcont_mean+Lcont_std,'b','LineWidth',1);
    %plot((1000*t1):1000*t2-1,Lcont_mean-Lcont_std,'b','LineWidth',1);
    
    plot((1000*t1):1000*t2-1,Hcont_mean,'r','LineWidth',3);
    %plot((1000*t1):1000*t2-1,Hcont_mean+Hcont_std,'r','LineWidth',1);
    %plot((1000*t1):1000*t2-1,Hcont_mean-Hcont_std,'r','LineWidth',1);
    
    title('Neural Activity','FontSize',16);
    xlabel(sprintf('Time from %s (ms)',alignopts{aligntype-5}),'FontSize',14);
    ylabel('Firing Rate (spikes/sec)','FontSize',14);
    
    %Lcount_count_m = Lcont_mean;
    
    OUTPUTS.fano.xc = (1000*t1):1000*t2-1;
    OUTPUTS.fano.Lc = Lcont_std.^2./Lcont_mean;
    OUTPUTS.fano.Hc = Hcont_std.^2./Hcont_mean;
    
end
if WindowBin > 0
    figure; hold on;

    Lwin = zeros(LL,length(winedges)-1);
    Hwin = zeros(LH,length(winedges)-1);
    for i = 1:length(winedges)-1
        Lwin(:,i) = sum(Lraster_ext(:,winedges(i):winedges(i+1)-1),2);
        Hwin(:,i) = sum(Hraster_ext(:,winedges(i):winedges(i+1)-1),2);
    end
    
    Lcount_m = mean(Lwin,1);
    Hcount_m = mean(Hwin,1);
    
    Lcount_std = std(Lwin,1);
    Hcount_std = std(Hwin,1);
    
    L_ff = (Lcount_std.^2)./Lcount_m;
    H_ff = (Hcount_std.^2)./Hcount_m;
    
    Lff_interp = interp1(wincents,L_ff,1:round(1000*(t2-t1))+2000);
    Hff_interp = interp1(wincents,H_ff,1:round(1000*(t2-t1))+2000);
    
    Lff_interp(1:1000) = []; Lff_interp(end-999:end) = [];
    Hff_interp(1:1000) = []; Hff_interp(end-999:end) = [];
    
    OUTPUTS.fano.xw = (1000*t1):1000*t2-1;
    OUTPUTS.fano.Lw = Lff_interp;
    OUTPUTS.fano.Hw = Hff_interp;
    
    %%% PLOTTING %%%
    LwinFR = Lcount_m.*1000./WindowBin;
    HwinFR = Hcount_m.*1000./WindowBin;
    Lwin_interp = interp1(wincents,LwinFR,1:round(1000*(t2-t1))+2000);
    Lwin_interp(1:1000) = []; Lwin_interp(end-999:end) = [];
    
    Hwin_interp = interp1(wincents,HwinFR,1:round(1000*(t2-t1))+2000);
    Hwin_interp(1:1000) = []; Hwin_interp(end-999:end) = [];
    
    % Do Plots of FR
    plot((1000*t1):1000*t2-1,Lwin_interp,'b','LineWidth',3);
    plot((1000*t1):1000*t2-1,Hwin_interp,'r','LineWidth',3);
    
    title('Neural Activity','FontSize',16);
    xlabel(sprintf('Time from %s (ms)',alignopts{aligntype-5}),'FontSize',14);
    ylabel('Firing Rate (spikes/sec)','FontSize',14);
      
end
    
%% Plot endpoint-based FR distribution

   
[ordered,Linds] = sort(Lpos_spikes,1);
Lpos_spikes_ordered = Lpos_spikes(Linds(:,1),:);
[ordered,Hinds] = sort(Hpos_spikes,1);
Hpos_spikes_ordered = Hpos_spikes(Hinds(:,1),:);

[Ln,Lbin] = histc(Lpos_spikes_ordered(:,1),xs);
[Hn,Hbin] = histc(Hpos_spikes_ordered(:,1),xs);

Lsumpos = zeros(length(xs)-1,1);
Lstdpos = zeros(length(xs)-1,1);
uniqueL = unique(Lbin);
for i = 1:length(uniqueL) 
    q = uniqueL(i);
    posinds = find(Lbin==q);
    Lsumpos(q) = nanmean(Lpos_spikes_ordered(posinds,2)); 
    Lstdpos(q) = nanstd(Lpos_spikes_ordered(posinds,2));
end
Hsumpos = zeros(length(xs)-1,1);
Hstdpos = zeros(length(xs)-1,1);
uniqueH = unique(Hbin);
for i = 1:length(uniqueH)  
    q = uniqueH(i);
    posinds = find(Hbin==q);
    Hsumpos(q) = nanmean(Hpos_spikes_ordered(posinds,2)); 
    Hstdpos(q) = nanstd(Hpos_spikes_ordered(posinds,2));       
end

Lsumpos(isnan(Lsumpos))=0;
Lstdpos(isnan(Lstdpos))=0;
Hsumpos(isnan(Hsumpos))=0;
Hstdpos(isnan(Hstdpos))=0;

plotx = zeros(1,2*length(xs)-2);
plotx(1:2:end) = xs(1:end-1);
plotx(2:2:end) = xs(2:end);

plotLsum = zeros(1,length(plotx));
plotHsum = zeros(1,length(plotx));

plotLsum(1:2:end) = Lsumpos';
plotLsum(2:2:end) = Lsumpos';
plotHsum(1:2:end) = Hsumpos';
plotHsum(2:2:end) = Hsumpos';    
    
%     figure; plot(plotx,plotLsum,'b');
%     hold on; plot(plotx,plotHsum,'r');
if Endpoint_rank == 1  && Errorbarbinplot == 1  
    figure; errorbar(.5*(xs(1:end-1)+xs(2:end)),Lsumpos',Lstdpos,'b');
    hold on; errorbar(.5*(xs(1:end-1)+xs(2:end)),Hsumpos',Hstdpos,'r');
        
end

%% Other useful metrics
if BlobPlot == 1
    low_isis = vertcat(L_ISI{:,:});
    high_isis = vertcat(H_ISI{:,:});

    L_nospikes = length(find(Lpos_spikes(:,2)==0));
    H_nospikes = length(find(Hpos_spikes(:,2)==0));

    fprintf('No Spikes: (%d %d)\n',L_nospikes,H_nospikes);

    range = (1000*t1):1000*t2-1;
    L_LP(L_LP(:,3)==0,3)=0.1;
    H_LP(H_LP(:,3)==0,3)=0.1;

    figure; subplot(2,1,1); hold on;
    for i = 1:size(L_LP,1)
        plot(L_LP(i,1),range(L_LP(i,2)),'b.','MarkerSize',L_LP(i,3));
        if Lpos_spikes(i,2)==0
            plot(repmat(L_LP(i,1),2,1),[0 10000],'b');   
        end
    end
    axis([xs(1) xs(end) 1000*t1-100 1000*t2+100]);
    subplot(2,1,2); hold on;
    for i = 1:size(H_LP,1)
        plot(H_LP(i,1),range(H_LP(i,2)),'r.','MarkerSize',H_LP(i,3));
        if Hpos_spikes(i,2)==0
            plot(repmat(H_LP(i,1),2,1),[0 10000],'r');   
        end

    end
    axis([xs(1) xs(end) 1000*t1-100 1000*t2+100]);
end

%% Binned movement PSTHs

if PosBinned == 1
    
    shift_contL = [Lpos_spikes(:,1) Lcont];
    shift_contH = [Hpos_spikes(:,1) Hcont];

    sc_L = shift_contL(Linds(:,1),:);
    sc_H = shift_contH(Hinds(:,1),:);

    SCLm = zeros(length(xs)-1,size(Lcont,2));
    SCLs = zeros(length(xs)-1,size(Lcont,2));
    for i = 1:length(uniqueL) 
        q = uniqueL(i);
        posinds = find(Lbin==q);
        SCLm(q,:) = nanmean(sc_L(posinds,2:end),1); 
        SCLs(q,:) = nanstd(sc_L(posinds,2:end),1);
    end
    SCHm = zeros(length(xs)-1,size(Hcont,2));
    SCHs = zeros(length(xs)-1,size(Hcont,2));
    for i = 1:length(uniqueH) 
        q = uniqueH(i);
        posinds = find(Hbin==q);
        SCHm(q,:) = nanmean(sc_H(posinds,2:end),1); 
        SCHs(q,:) = nanstd(sc_H(posinds,2:end),1);
    end
    
    figure; hold on;
    for i = 1:length(xs)-1
        plot((1000*t1):1000*t2-1,SCLm(i,:),'Linewidth',4,...
            'Color',[i/(length(xs)+1) i/(length(xs)+1) 1]);
        
        plot((1000*t1):1000*t2-1,SCHm(i,:),'Linewidth',4,...
            'Color',[1 i/(length(xs)+1) i/(length(xs)+1)]);
    end
        

end