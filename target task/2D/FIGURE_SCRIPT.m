%% FR Gain RAW
%good_units = [15.2 37.2 41.2 49.2 51.1 53.1 59.1 63.1 82.3 88.2]
good_indices = 104;%[25 74 84 104 108 111 125 136 166 175];
unit_list = good_indices;
for UNITIND = 1:length(unit_list)
    testunit = unit_list(UNITIND);

    coinds = find(abs(circ_dist(alldays(1).tt(:,10),pi/2)) < 0.3); 
    uL = find(abs(circ_dist(alldays(2).tt(:,10),pi/2)) < 0.1 & alldays(2).tt(:,3)==max(alldays(2).tt(:,3)));
    matchedL_int = randperm(length(uL));
    matchedL = matchedL_int(1:length(coinds));
    %matchedL = randperm(length(uL),length(coinds));
    uLinds = uL(matchedL);

    uH = find(abs(circ_dist(alldays(2).tt(:,10),pi/2)) < 0.1 & alldays(2).tt(:,3)==min(alldays(2).tt(:,3)));
    matchedH_int = randperm(length(uH));
    matchedH = matchedH_int(1:length(coinds));
    %matchedH = randperm(length(uH),length(coinds));
    uHinds = uH(matchedH);

    uM = find(abs(circ_dist(alldays(2).tt(:,10),pi/2)) < 0.1 & alldays(2).tt(:,3)==5);
    matchedM_int = randperm(length(uM));
    matchedM = matchedM_int(1:length(coinds));
    %matchedM = randperm(length(uM),length(coinds));
    uMinds = uM(matchedM);

    cott = alldays(1).tt(coinds,:);
    Ltt = alldays(2).tt(uLinds,:);
    Htt = alldays(2).tt(uHinds,:);
    Mtt = alldays(2).tt(uMinds,:);

    TT = {cott,Ltt,Htt};

    colorp = {'k','b','r'};
    figure; hold on;
    RAST = cell(3,2);
    PSTH = cell(3,2);
    for i = 1:3

        nT = size(TT{i},1);
        day_spikes = alldays(1).PMd_units;
        day_tt = TT{i}; %alldays(1).tt;

        t1_t2 = [-.2 .8]; % Time (s) before and after event
        rxT = repmat(t1_t2(1):0.001:(t1_t2(2)-0.001),nT,1);
        pxT = [-0.1:0.2:0.7];
        alignment = 'target'; % event to align with ('target','go','reward')
        [rast1,~,~] = raster_plot(day_spikes{testunit},day_tt,t1_t2,...
            alignment,0,'trials');

        rast = rast1{1} + length(cott)*(i-1);
        psth1 = bin_array(rast1{1}./rast1{1},1,5,'sum');

        RAST{i,1} = rast; PSTH{i,1} = (psth1/nT)/0.2;

        t1_t2 = [0 .6]; % Time (s) before and after event
        rxG = repmat(t1_t2(1):0.001:(t1_t2(2)-0.001),nT,1) + 1;
        pxG = [1.1:0.2:1.5];
        alignment = 'go'; % event to align with ('target','go','reward')
        [rast2,~,~] = raster_plot(day_spikes{testunit},day_tt,t1_t2,...
            alignment,0,'trials');

        rast = rast2{1} + length(cott)*(i-1);
        psth2 = bin_array(rast2{1}./rast2{1},1,3,'sum');

        RAST{i,2} = rast; PSTH{i,2} = (psth2/nT)/0.2;
    end

    PSTHMAX = max(horzcat(PSTH{:}))+1;
    for i = 1:3 

        plot(rxT(~isnan(RAST{i,1})),RAST{i,1}(~isnan(RAST{i,1}))./(3*nT./PSTHMAX),'.','Color','k'); 
        plot(rxG(~isnan(RAST{i,2})),RAST{i,2}(~isnan(RAST{i,2}))./(3*nT./PSTHMAX),'.','Color','k');
       
        plot(pxT,PSTH{i,1},'Color',colorp{i},'LineWidth',4); plot(pxG,PSTH{i,2},'Color',colorp{i},'LineWidth',4);

        patch([-0.2 1.6 1.6 -0.2],((i-1)*nT + nT.*[0.01 0.01 1.01 1.01])./(3*nT./PSTHMAX),colorp{i},'FaceAlpha',1,'EdgeAlpha',1);
    end
    plot([0 0],PSTHMAX.*[0,1.01],'k--');
    plot([0.8 0.8],PSTHMAX.*[0,1.01],'k--');
    plot([1 1],PSTHMAX.*[0,1.01],'k--');
    xlabel('Time from Target ON (s)','FontSize',14);
    ylabel('Firing Rate (spikes/sec)','FontSize',14);
    title(sprintf('Unit %.1f',alldays(1).PMd_units{testunit}(1)),'FontSize',16);
    ylim([0 (2*nT + nT*1.01)./(3*nT./PSTHMAX)]);

    pause; close;
end

%% Tuning Width RAW
spreaduns = [42 73 128];
figure; hold on;
for unis = 1:3
    testunit = spreaduns(unis);
    subplot(1,3,unis); hold on;
    coinds = find(abs(circ_dist(alldays(1).tt(:,10),pi/2)) < pi/5); 
    uL = find(abs(circ_dist(alldays(2).tt(:,10),pi/2)) < 0.1 & alldays(2).tt(:,3)==max(alldays(2).tt(:,3)));
    matchedL_int = randperm(length(uL));
    matchedL = matchedL_int(1:length(coinds));
    uLinds = uL(matchedL);

    uH = find(abs(circ_dist(alldays(2).tt(:,10),pi/2)) < 0.1 & alldays(2).tt(:,3)==min(alldays(2).tt(:,3)));
    matchedH_int = randperm(length(uH));
    matchedH = matchedH_int(1:length(coinds));
    uHinds = uH(matchedH);

    uM = find(abs(circ_dist(alldays(2).tt(:,10),pi/2)) < 0.1 & alldays(2).tt(:,3)==5);
    matchedM_int = randperm(length(uM));
    matchedM = matchedM_int(1:length(coinds));
    uMinds = uM(matchedM);

    cott = alldays(1).tt(coinds,:);
    Ltt = alldays(2).tt(uLinds,:);
    Htt = alldays(2).tt(uHinds,:);
    Mtt = alldays(2).tt(uMinds,:);

    TT = {cott,Ltt,Htt};

    colorp = {'k','b','r'};

    RAST = cell(3,2);
    PSTH = cell(3,2);
    for i = 1:3

        nT = size(TT{i},1);
        day_spikes = alldays(1).PMd_units;
        day_tt = TT{i}; %alldays(1).tt;

        t1_t2 = [-.2 .8]; % Time (s) before and after event
        rxT = repmat(t1_t2(1):0.001:(t1_t2(2)-0.001),nT,1);
        pxT = [-0.1:0.2:0.7];
        alignment = 'target'; % event to align with ('target','go','reward')
        [rast1,~,~] = raster_plot(day_spikes{testunit},day_tt,t1_t2,...
            alignment,0,'trials');

        rast = rast1{1} + length(cott)*(i-1);
        psth1 = bin_array(rast1{1}./rast1{1},1,5,'sum');

        RAST{i,1} = rast; PSTH{i,1} = (psth1/nT)/0.2;

        t1_t2 = [0 .6]; % Time (s) before and after event
        rxG = repmat(t1_t2(1):0.001:(t1_t2(2)-0.001),nT,1) + 1;
        pxG = [1.1:0.2:1.5];
        alignment = 'go'; % event to align with ('target','go','reward')
        [rast2,~,~] = raster_plot(day_spikes{testunit},day_tt,t1_t2,...
            alignment,0,'trials');

        rast = rast2{1} + length(cott)*(i-1);
        psth2 = bin_array(rast2{1}./rast2{1},1,3,'sum');

        RAST{i,2} = rast; PSTH{i,2} = (psth2/nT)/0.2;
    end

    PSTHMAX = max(horzcat(PSTH{:}))+1;
    for i = 1:3 

        plot(rxT,RAST{i,1}./(3*nT./PSTHMAX),'.','Color','k'); plot(rxG,RAST{i,2}./(3*nT./PSTHMAX),'.','Color','k');
        plot(pxT,PSTH{i,1},'Color',colorp{i},'LineWidth',4); plot(pxG,PSTH{i,2},'Color',colorp{i},'LineWidth',4);

        patch([-0.2 1.6 1.6 -0.2],((i-1)*nT + nT.*[0.01 0.01 1.01 1.01])./(3*nT./PSTHMAX),colorp{i},'FaceAlpha',0.25,'EdgeAlpha',0);
    end
    plot([0 0],PSTHMAX.*[0,1.01],'k--');
    plot([0.8 0.8],PSTHMAX.*[0,1.01],'k--');
    plot([1 1],PSTHMAX.*[0,1.01],'k--');
    xlabel('Time from Target ON (s)','FontSize',14);
    ylabel('Firing Rate (spikes/sec)','FontSize',14);
    title(sprintf('Unit %.1f',PMd_units{testunit}(1)),'FontSize',16);
    ylim([0 (2*nT + nT*1.01)./(3*nT./PSTHMAX)]);
end

%%

neuron_list =[111 166 90 25 74 84 104 108  125 136  175];
for i = 1:length(neuron_list)
    
    unis = neuron_list(i);
    
    t1_t2 = [200 800];
    
    co_units = alldays(1).PMd_units;
    
    [cents_co,tune_co,~,~,ind_co] = co_tuning(co_units,alldays(1).tt(:,10),...
                alldays(1).tt,t1_t2(1),t1_t2(2),'target',unis);
            
    test_units = alldays(2).PMd_units;
    
    [cents_teL,tune_teL,~,~,ind_teL] = co_tuning(test_units,Ltt(:,10),...
                Ltt,t1_t2(1),t1_t2(2),'target',unis);
    [cents_teM,tune_teM,~,~,ind_teM] = co_tuning(test_units,Mtt(:,10),...
               Mtt,t1_t2(1),t1_t2(2),'target',unis);
    [cents_teH,tune_teH,~,~,ind_teH] = co_tuning(test_units,Htt(:,10),...
                Htt,t1_t2(1),t1_t2(2),'target',unis);
    
    triallengthsL = cellfun(@length,ind_teL);
    bad_dirsL = find(triallengths < 20);
        
    triallengthsM = cellfun(@length,ind_teM);
    bad_dirsM = find(triallengths < 20);
        
    triallengthsH = cellfun(@length,ind_teH);
    bad_dirsH = find(triallengths < 20);
            
    tune_co(tune_co==0) = NaN;
    tune_teL(tune_teL==0) = NaN;
    tune_teL(bad_dirsL) = NaN;
    
    tune_teM(tune_teM==0) = NaN;
    tune_teM(bad_dirsM) = NaN;
    
    tune_teH(tune_teH==0) = NaN;
    tune_teH(bad_dirsH) = NaN;
            
    polygon_tuning(cents_co,tune_co,1);
    polygon_tuning(cents_teL,tune_teL,2);
    polygon_tuning(cents_teM,tune_teM,3);
    polygon_tuning(cents_teH,tune_teH,4);
    
    pause; close;
end

%% Fdiff with tuning curves
%lowinds = [55 65 121 135 158 173 174];
%highinds = [25 100 108 111];
spreaduns = [lowinds(3) highinds(4)];
figure; hold on;
for unis = 1:length(spreaduns)
    testunit = spreaduns(unis);
    subplot(2,length(spreaduns),unis); hold on;
    coinds = find(abs(circ_dist(alldays(1).tt(:,10),pi/2)) < pi/5); 
    uL = find(abs(circ_dist(alldays(2).tt(:,10),pi/2)) < 0.1 & alldays(2).tt(:,3)==max(alldays(2).tt(:,3)));
    matchedL_int = randperm(length(uL));
    matchedL = matchedL_int(1:length(coinds));
    uLinds = uL;%(matchedL);

    uH = find(abs(circ_dist(alldays(2).tt(:,10),pi/2)) < 0.1 & alldays(2).tt(:,3)==min(alldays(2).tt(:,3)));
    matchedH_int = randperm(length(uH));
    matchedH = matchedH_int(1:length(coinds));
    uHinds = uH;%(matchedH);

    uM = find(abs(circ_dist(alldays(2).tt(:,10),pi/2)) < 0.1 & alldays(2).tt(:,3)==5);
    matchedM_int = randperm(length(uM));
    matchedM = matchedM_int(1:length(coinds));
    uMinds = uM;%(matchedM);

    cott = alldays(1).tt(coinds,:);
    Ltt = alldays(2).tt(uLinds,:);
    Htt = alldays(2).tt(uHinds,:);
    Mtt = alldays(2).tt(uMinds,:);

    TT = {cott,Ltt,Htt};

    colorp = {'k','b','r'};

    RAST = cell(3,2);
    [PSTH,PSTH_L,PSTH_H] = deal(cell(3,2));
    for i = 1:3

        nT = size(TT{i},1);
        day_spikes = alldays(1).PMd_units;
        day_tt = TT{i}; %alldays(1).tt;

        t1_t2 = [-.2 .8]; % Time (s) before and after event
        rxT = repmat(t1_t2(1):0.001:(t1_t2(2)-0.001),nT,1);
        pxT = [-0.1:0.2:0.7];
        alignment = 'target'; % event to align with ('target','go','reward')
        [rast1,~,~] = raster_plot(day_spikes{testunit},day_tt,t1_t2,...
            alignment,0,'trials');

        rast = rast1{1} + length(cott)*(i-1);
        psth1 = bin_array(rast1{1}./rast1{1},1,5,'sum');
        
        [psth1L,psth1H] = boot_bounds(1000,@(x) bin_array(x./x,1,5,'sum'), rast1{1}, 2.5,97.5);

        RAST{i,1} = rast; PSTH{i,1} = (psth1/nT)/0.2;
        PSTH_L{i,1} = (psth1L'/nT)/0.2;
        PSTH_H{i,1} = (psth1H'/nT)/0.2;

        t1_t2 = [0 .6]; % Time (s) before and after event
        rxG = repmat(t1_t2(1):0.001:(t1_t2(2)-0.001),nT,1) + 1;
        pxG = [1.1:0.2:1.5];
        alignment = 'go'; % event to align with ('target','go','reward')
        [rast2,~,~] = raster_plot(day_spikes{testunit},day_tt,t1_t2,...
            alignment,0,'trials');

        rast = rast2{1} + length(cott)*(i-1);
        psth2 = bin_array(rast2{1}./rast2{1},1,3,'sum');
        
        [psth2L,psth2H] = boot_bounds(1000,@(x) bin_array(x./x,1,3,'sum'), rast2{1}, 2.5,97.5);

        RAST{i,2} = rast; PSTH{i,2} = (psth2/nT)/0.2;
        PSTH_L{i,2} = (psth2L'/nT)/0.2;
        PSTH_H{i,2} = (psth2H'/nT)/0.2;
    end

    PSTHMAX = max(horzcat(PSTH_H{:}))+1;
    for i = 1:3 

        %plot(rxT,RAST{i,1}./(3*nT./PSTHMAX),'.','Color','k'); plot(rxG,RAST{i,2}./(3*nT./PSTHMAX),'.','Color','k');
        plot(pxT,PSTH{i,1},'Color',colorp{i},'LineWidth',3); plot(pxG,PSTH{i,2},'Color',colorp{i},'LineWidth',3);
        patch([pxT fliplr(pxT)],[PSTH_L{i,1} fliplr(PSTH_H{i,1})],colorp{i},'FaceAlpha',0.25,'EdgeAlpha',0);
        patch([pxG fliplr(pxG)],[PSTH_L{i,2} fliplr(PSTH_H{i,2})],colorp{i},'FaceAlpha',0.25,'EdgeAlpha',0);

        %patch([-0.2 1.6 1.6 -0.2],((i-1)*nT + nT.*[0.01 0.01 1.01 1.01])./(3*nT./PSTHMAX),colorp{i},'FaceAlpha',0.25,'EdgeAlpha',0);
    end
    plot([0 0],PSTHMAX.*[0,1.01],'k--');
    plot([0.8 0.8],PSTHMAX.*[0,1.01],'k--');
    plot([1 1],PSTHMAX.*[0,1.01],'k--');
    xlabel('Time from Target ON (s)','FontSize',14);
    ylabel('Firing Rate (spikes/sec)','FontSize',14);
    title(sprintf('Unit %.1f',PMd_units{testunit}(1)),'FontSize',16);
    ylim([0 PSTHMAX]);
    xlim([-0.2 1.6]);
    
    
    subplot(2,length(spreaduns),unis+length(spreaduns)); hold on;
    TC = neurons{1}{testunit}.tuning;
    
    move_region = find(abs(circ_dist(wrapped_cents,pi/2))<0.1);
    
    plot(wrapped_cents,TC,'k','LineWidth',2);
    patch([pi/2-0.1 pi/2+0.1 (fliplr(wrapped_cents(move_region)))],...
        [0.89*min(TC) 0.89*min(TC) (flipud(TC(move_region))')],'k','FaceAlpha',0.5);
    
%     plot(wrapped_cents(abs(circ_dist(wrapped_cents,pi/2))<0.1),...
%         TC(abs(circ_dist(wrapped_cents,pi/2))<0.1),'k','LineWidth',5);
    
    ylim([0.9*min(TC) 1.1*max(TC)]);
    xlim([0 2*pi]);
    
    xlabel('\theta','FontSize',14);
    ylabel('Firing Rate (/sec)','FontSize',14);
end

%% Plain old rasters
%good_units = [15.2 37.2 41.2 49.2 51.1 53.1 59.1 63.1 82.3 88.2]
%good_indices = [25 74 84 104 108 111 125 136 166 175];
unit_list = 19;%1:length(PMd_ids);
for UNITIND = 1:length(unit_list)
    testunit = unit_list(UNITIND);

    Htt = alldays(2).tt;
    TT = Htt; TT(:,3)=1;
    
    figure; hold on;
    RAST = cell(1,2);
    PSTH = cell(1,2);
    i = 1;

    nT = size(TT,1);
    day_spikes = alldays(1).PMd_units;
    day_tt = TT; %alldays(1).tt;

    t1_t2 = [-1 .8]; % Time (s) before and after event
    rxT = repmat(t1_t2(1):0.001:(t1_t2(2)-0.001),nT,1);
    pxT = t1_t2(1):0.2:t1_t2(2);
    dT = diff(pxT(1:2));
    alignment = 'target'; % event to align with ('target','go','reward')
    [rast1,~,~] = raster_plot(day_spikes{testunit},day_tt,t1_t2,...
        alignment,0,'trials');

    rast = rast1{1};
    psth1 = bin_array(rast./rast,1,5,'sum');

    RAST{1} = rast; PSTH{1} = (psth1/nT)/dT;

    t1_t2 = [0 1.5]; % Time (s) before and after event
    rxG = repmat(t1_t2(1):0.001:(t1_t2(2)-0.001),nT,1) + 1;
    pxG = t1_t2(1):0.2:t1_t2(2);
    alignment = 'go'; % event to align with ('target','go','reward')
    [rast2,~,~] = raster_plot(day_spikes{testunit},day_tt,t1_t2,...
        alignment,0,'trials');

    rast = rast2{1};
    psth2 = bin_array(rast./rast,1,3,'sum');

    RAST{2} = rast; PSTH{2} = (psth2/nT)/dT;
    
    PSTHMAX = size(RAST{1},1);

    %plot(rxT,RAST{1},'.','Color','k','MarkerSize',10); plot(rxG,RAST{2},'.','Color','k','MarkerSize',10);
    plot(rxT(~isnan(RAST{1})),RAST{1}(~isnan(RAST{1})),'.','Color','k','MarkerSize',10); 
    plot(rxG(~isnan(RAST{2})),RAST{2}(~isnan(RAST{2})),'.','Color','k','MarkerSize',10);
    %plot(pxT,PSTH{1},'Color','k','LineWidth',4); plot(pxG,PSTH{2},'Color','k','LineWidth',4);
    %patch([-0.2 1.6 1.6 -0.2],((i-1)*nT + nT.*[0.01 0.01 1.01 1.01])./(3*nT./PSTHMAX),colorp{i},'FaceAlpha',0.25,'EdgeAlpha',0);
    
    plot([0 0],PSTHMAX.*[0,1.01],'b--','LineWidth',3);
    plot([0.8 0.8],PSTHMAX.*[0,1.01],'b--','LineWidth',3);
    plot([1 1],PSTHMAX.*[0,1.01],'b--','LineWidth',3);
    xlabel('Time from Target ON (s)','FontSize',14);
    ylabel('Trial','FontSize',14);
    title(sprintf('Unit %.1f',PMd_units{testunit}(1)),'FontSize',16);
    ylim([0 PSTHMAX]);

end

%% Distributions

thets = -pi:0.01:pi;
pridist = circ_vmpdf(thets,pi/2,25);
likdist1 = circ_vmpdf(thets,0,1);
likdist2 = circ_vmpdf(thets,0,5);
likdist3 = circ_vmpdf(thets,0,50);

figure;
plot(thets,pridist,'k','LineWidth',4); ylim([0 3]); xlim([-pi pi]); box off;

figure; plot(thets,likdist1,'r','LineWidth',4); ylim([0 3]); xlim([-pi pi]); box off;
figure; plot(thets,likdist2,'g','LineWidth',4); ylim([0 3]); xlim([-pi pi]); box off;
figure; plot(thets,likdist3,'b','LineWidth',4); ylim([0 3]); xlim([-pi pi]); box off;

%% Speeds
Linds = find(alldays(2).tt(:,3)==100);
%Minds = find(alldays(2).tt(:,3)==5);
Hinds = find(alldays(2).tt(:,3)==1);

Ltt = alldays(2).tt(Linds,:); [speedsL_t,xL_t,yL_t] = kin_exam(alldays(2).bdfM,Ltt,800,5);
%Mtt = alldays(2).tt(Minds,:); [speedsM_t,xM_t,yM_t] = kin_exam(alldays(2).bdfM,Mtt,800,5);
Htt = alldays(2).tt(Hinds,:); [speedsH_t,xH_t,yH_t] = kin_exam(alldays(2).bdfM,Htt,800,5);

[speedsL_g,xL_g,yL_g] = kin_exam(alldays(2).bdfM,Ltt,1000,6);
%[speedsM_g,xM_g,yM_g] = kin_exam(alldays(2).bdfM,Mtt,1000,6);
[speedsH_g,xH_g,yH_g] = kin_exam(alldays(2).bdfM,Htt,1000,6);

figure; hold on;
plot([0:800],mean(speedsL_t,1),'b','LineWidth',2);
%plot([0:800],mean(speedsM_t,1),'g','LineWidth',2);
plot([0:800],mean(speedsH_t,1),'r','LineWidth',2);

[LLT, ULT] = boot_bounds(1000,@mean,speedsL_t,2.5,97.5);
%[LMT, UMT] = boot_bounds(1000,@mean,speedsM_t,2.5,97.5);
[LHT, UHT] = boot_bounds(1000,@mean,speedsH_t,2.5,97.5);

patch([0:800 800:-1:0],[LLT' fliplr(ULT')],'b','FaceAlpha',0.25,'EdgeAlpha',0);
%patch([0:800 800:-1:0],[LMT' fliplr(UMT')],'g','FaceAlpha',0.25,'EdgeAlpha',0);
patch([0:800 800:-1:0],[LHT' fliplr(UHT')],'r','FaceAlpha',0.25,'EdgeAlpha',0);

plot([1000:2000],mean(speedsL_g,1),'b','LineWidth',2);
%plot([1000:2000],mean(speedsM_g,1),'g','LineWidth',2);
plot([1000:2000],mean(speedsH_g,1),'r','LineWidth',2);

[LLG, ULG] = boot_bounds(1000,@mean,speedsL_g,2.5,97.5);
%[LMG, UMG] = boot_bounds(1000,@mean,speedsM_g,2.5,97.5);
[LHG, UHG] = boot_bounds(1000,@mean,speedsH_g,2.5,97.5);

patch([1000:2000 2000:-1:1000],[LLG' fliplr(ULG')],'b','FaceAlpha',0.25,'EdgeAlpha',0);
%patch([1000:2000 2000:-1:1000],[LMG' fliplr(UMG')],'g','FaceAlpha',0.25,'EdgeAlpha',0);
patch([1000:2000 2000:-1:1000],[LHG' fliplr(UHG')],'r','FaceAlpha',0.25,'EdgeAlpha',0);

plot([800 1000; 800 1000],[0 0; 25 25],'k--');

%% FR Gain RAW
%good_units = [15.2 37.2 41.2 49.2 51.1 53.1 59.1 63.1 82.3 88.2]
good_indices = 104;%[25 74 84 104 108 111 125 136 166 175];
unit_list = good_indices;
for UNITIND = 1:length(unit_list)
    testunit = unit_list(UNITIND);

    coinds = find(abs(circ_dist(alldays(1).tt(:,10),pi/2)) < 0.3); 
    uL = find(abs(circ_dist(alldays(2).tt(:,10),pi/2)) < 0.1 & alldays(2).tt(:,3)==max(alldays(2).tt(:,3)));
    matchedL_int = randperm(length(uL));
    matchedL = matchedL_int(1:length(coinds));
    %matchedL = randperm(length(uL),length(coinds));
    uLinds = uL(matchedL);

    uH = find(abs(circ_dist(alldays(2).tt(:,10),pi/2)) < 0.1 & alldays(2).tt(:,3)==min(alldays(2).tt(:,3)));
    matchedH_int = randperm(length(uH));
    matchedH = matchedH_int(1:length(coinds));
    %matchedH = randperm(length(uH),length(coinds));
    uHinds = uH(matchedH);

    uM = find(abs(circ_dist(alldays(2).tt(:,10),pi/2)) < 0.1 & alldays(2).tt(:,3)==5);
    matchedM_int = randperm(length(uM));
    matchedM = matchedM_int(1:length(coinds));
    %matchedM = randperm(length(uM),length(coinds));
    uMinds = uM(matchedM);

    cott = alldays(1).tt(coinds,:);
    Ltt = alldays(2).tt(uLinds,:);
    Htt = alldays(2).tt(uHinds,:);
    Mtt = alldays(2).tt(uMinds,:);

    TT = {cott,Ltt,Htt};

    colorp = {'k','b','r'};
    figure; hold on;
    [PSTH,PSTH_L,PSTH_U] = deal(cell(3,2));
    for i = 1:3

        nT = size(TT{i},1);
        day_spikes = alldays(1).PMd_units;
        day_tt = TT{i}; %alldays(1).tt;

        t1_t2 = [-.2 .8]; % Time (s) before and after event
        pxT = [-0.1:0.2:0.7];
        alignment = 'target'; % event to align with ('target','go','reward')
        [rast1,~,~] = raster_plot(day_spikes{testunit},day_tt,t1_t2,...
            alignment,0,'trials');

        psth1 = bin_array(rast1{1}./rast1{1},1,5,'sum');
        [psth1_l, psth1_u] = boot_bounds(1000,@(x) bin_array(x,1,5,'sum'),rast1{1}./rast1{1},2.5,97.5);

        PSTH{i,1} = (psth1/nT)/0.2;
        PSTH_L{i,1} = (psth1_l/nT)'/0.2;
        PSTH_U{i,1} = (psth1_u/nT)'/0.2;


        t1_t2 = [0 .6]; % Time (s) before and after event
        pxG = [1.1:0.2:1.5];
        alignment = 'go'; % event to align with ('target','go','reward')
        [rast2,~,~] = raster_plot(day_spikes{testunit},day_tt,t1_t2,...
            alignment,0,'trials');

        psth2 = bin_array(rast2{1}./rast2{1},1,3,'sum');
        [psth2_l, psth2_u] = boot_bounds(1000,@(x) bin_array(x,1,3,'sum'),rast2{1}./rast2{1},2.5,97.5);

        PSTH{i,2} = (psth2/nT)/0.2;
        PSTH_L{i,2} = (psth2_l/nT)'/0.2;
        PSTH_U{i,2} = (psth2_u/nT)'/0.2;
    end

    PSTHMAX = max(horzcat(PSTH{:}))+1;
    for i = 1:3 

        plot(pxT,PSTH{i,1},'Color',colorp{i},'LineWidth',4); 
        plot(pxG,PSTH{i,2},'Color',colorp{i},'LineWidth',4);
        
        patch([pxT fliplr(pxT)],[PSTH_L{i,1} fliplr(PSTH_U{i,1})],colorp{i},'FaceAlpha',0.25,'EdgeAlpha',0);
        patch([pxG fliplr(pxG)],[PSTH_L{i,2} fliplr(PSTH_U{i,2})],colorp{i},'FaceAlpha',0.25,'EdgeAlpha',0);

    end
    plot([0 0],PSTHMAX.*[0,1.01],'k--');
    plot([0.8 0.8],PSTHMAX.*[0,1.01],'k--');
    plot([1 1],PSTHMAX.*[0,1.01],'k--');
    xlabel('Time from Target ON (s)','FontSize',14);
    ylabel('Firing Rate (spikes/sec)','FontSize',14);
    title(sprintf('Unit %.1f',alldays(1).PMd_units{testunit}(1)),'FontSize',16);
    ylim([0 (2*nT + nT*1.01)./(3*nT./PSTHMAX)]);

end


            