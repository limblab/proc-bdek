%% RT condition
meanfunc = @(x) [nanmean(x(x(:,2)==min(x(:,2)),1)),nanmean(x(x(:,2)==max(x(:,2)),1))];
medianfunc = @(x) [nanmedian(x(x(:,2)==min(x(:,2)),1)),nanmedian(x(x(:,2)==max(x(:,2)),1))];
[LIopp,ISopp,LIsame,ISsame] = deal(cell(length(LI),1));
for i = 1:length(LI)
    
    mrts = medianfunc(REACTTIME{i});
    signdiff = sign(diff(mrts));
    
    if signdiff == -1
        LIopp{i}{1} = LI{i}{1}(REACTTIME{i}(LI{i}{1}) > mrts(1));
        LIopp{i}{2} = LI{i}{2}(REACTTIME{i}(LI{i}{2}) < mrts(2) & REACTTIME{i}(LI{i}{2}) > 0);
        
        LIsame{i}{1} = LI{i}{1}(REACTTIME{i}(LI{i}{1}) < mrts(1) & REACTTIME{i}(LI{i}{1}) > 0);
        LIsame{i}{2} = LI{i}{2}(REACTTIME{i}(LI{i}{2}) > mrts(2));
    else
        LIopp{i}{1} = LI{i}{1}(REACTTIME{i}(LI{i}{1}) < mrts(1) & REACTTIME{i}(LI{i}{1}) > 0);
        LIopp{i}{2} = LI{i}{2}(REACTTIME{i}(LI{i}{2}) > mrts(2));
        
        LIsame{i}{1} = LI{i}{1}(REACTTIME{i}(LI{i}{1}) > mrts(1));
        LIsame{i}{2} = LI{i}{2}(REACTTIME{i}(LI{i}{2}) < mrts(2) & REACTTIME{i}(LI{i}{2}) > 0);
    end
    
    ISopp{i} = sortrows(vertcat(LIopp{i}{:}));
    ISsame{i} = sortrows(vertcat(LIsame{i}{:}));
end

%% flip peakspeed
meanfunc = @(x) [nanmean(x(x(:,2)==min(x(:,2)),1)),nanmean(x(x(:,2)==max(x(:,2)),1))];
medianfunc = @(x) [nanmedian(x(x(:,2)==min(x(:,2)),1)),nanmedian(x(x(:,2)==max(x(:,2)),1))];
[LIopp,ISopp,LIsame,ISsame] = deal(cell(length(LI),1));
for i = 1:length(LI)
    
    mrts = medianfunc(PEAKSPEED{i});
    signdiff = sign(diff(mrts));
    
    if signdiff == -1
        LIopp{i}{1} = LI{i}{1}(PEAKSPEED{i}(LI{i}{1}) > mrts(1));
        LIopp{i}{2} = LI{i}{2}(PEAKSPEED{i}(LI{i}{2}) < mrts(2) & PEAKSPEED{i}(LI{i}{2}) > 0);
        
        LIsame{i}{1} = LI{i}{1}(PEAKSPEED{i}(LI{i}{1}) < mrts(1) & PEAKSPEED{i}(LI{i}{1}) > 0);
        LIsame{i}{2} = LI{i}{2}(PEAKSPEED{i}(LI{i}{2}) > mrts(2));
    else
        LIopp{i}{1} = LI{i}{1}(PEAKSPEED{i}(LI{i}{1}) < mrts(1) & PEAKSPEED{i}(LI{i}{1}) > 0);
        LIopp{i}{2} = LI{i}{2}(PEAKSPEED{i}(LI{i}{2}) > mrts(2));
        
        LIsame{i}{1} = LI{i}{1}(PEAKSPEED{i}(LI{i}{1}) > mrts(1));
        LIsame{i}{2} = LI{i}{2}(PEAKSPEED{i}(LI{i}{2}) < mrts(2) & PEAKSPEED{i}(LI{i}{2}) > 0);
    end
    
    ISopp{i} = sortrows(vertcat(LIopp{i}{:}));
    ISsame{i} = sortrows(vertcat(LIsame{i}{:}));
end

%% RT all
dRTS = nan(length(REACTTIME),1);
dRTS_bnd = nan(length(REACTTIME),2);

dmeanfunc = @(x) nanmean(x(x(:,2)==min(x(:,2)),1))-nanmean(x(x(:,2)==max(x(:,2)),1));
Indall = cell(length(REACTTIME),1);
for i = 1:length(REACTTIME)
    
    Indall{i} = (REACTTIME{i}(:,1) < nanmedian(REACTTIME{i}(:,1)))*2 + ...
                    (REACTTIME{i}(:,1) > nanmedian(REACTTIME{i}(:,1)))*1;
                
    dRTS(i) = dmeanfunc(REACTTIME{i});
    
    [dRTS_bnd(i,1), dRTS_bnd(i,2)] = boot_bounds(1000,dmeanfunc,REACTTIME{i},2.5,97.5);

end


%% RT
figure; hold on; 
[sameRT] = deal(zeros(length(REACTTIME),1));
diffRT = zeros(length(REACTTIME),2);
for i = 1:length(REACTTIME)
    
    lindsF = @(x) find(x(:,2)==max(x(:,2)));
    hindsF = @(x) find(x(:,2)==min(x(:,2)));
    
    linds = lindsF(REACTTIME{i});
    hinds = hindsF(REACTTIME{i});
    
    rtL = REACTTIME{i}(linds,1);
    rtH = REACTTIME{i}(hinds,1);
    
    [sameRT(i),~,diffRT(i,2:3)] = ttest2(rtH,rtL);
    diffRT(i,1) = nanmean(rtH)-nanmean(rtL);
    if sameRT(i)==0
        mrk = '.';
    else
        mrk = 'o';
    end
    plot(i,diffRT(i,1),mrk,'Color','k');
    plot([i i],diffRT(i,2:3),'k');
    
%     plot(i,nanmean(rtL),mrk,'Color','b'); plot(i+0.2,nanmean(rtH),mrk,'Color','r');
end

%% Reaction time Flipped
figure; hold on; 
[sameRTf,stdsub] = deal(zeros(length(REACTTIME),1));
[diffRTf] = deal(zeros(length(REACTTIME),3));
[RT,newLI,rtF] = deal(cell(length(REACTTIME),1));
for i = G
    
    [okinds,stdsub(i)] = flip_dists(REACTTIME{G(i)});
    RT{i} = REACTTIME{G(i)}(okinds,:);
    
    newLI{i}{1} = okinds(REACTTIME{G(i)}(okinds,2)==max(REACTTIME{G(i)}(:,2)));
    newLI{i}{2} = okinds(REACTTIME{G(i)}(okinds,2)==min(REACTTIME{G(i)}(:,2)));
    
    lindsF = @(x) find(x(:,2)==max(x(:,2)));
    hindsF = @(x) find(x(:,2)==min(x(:,2)));
    
    linds = lindsF(RT{i});
    hinds = hindsF(RT{i});
    
    rtF{i,1} = RT{i}(linds,1);
    rtF{i,2} = RT{i}(hinds,1);
    
    [sameRTf(i),~,diffRTf(i,2:3)] = ttest2(rtF{i,2},rtF{i,1});
    diffRTf(i,1) = nanmean(rtF{i,2})-nanmean(rtF{i,1});
    if sameRTf(i)==0
        mrk = '.';
    else
        mrk = 'o';
    end
    plot(i,diffRTf(i,1),mrk,'Color','k');
    plot([i i],diffRTf(i,2:3),'k');
%     plot(i,nanmean(psL),mrk,'Color','b'); plot(i+0.2,nanmean(psH),mrk,'Color','r');
end

%% Peak speed
figure; hold on; 
[samePS] = deal(zeros(length(G),1));
diffPS = zeros(length(G),3);
for i = G%1:length(PEAKSPEED)
    
    lindsF = @(x) find(x(:,2)==max(x(:,2)));
    hindsF = @(x) find(x(:,2)==min(x(:,2)));
    
    linds = lindsF(PEAKSPEED{G(i)});
    hinds = hindsF(PEAKSPEED{G(i)});
    
    
    psL = PEAKSPEED{G(i)}(linds,1);
    psH = PEAKSPEED{G(i)}(hinds,1);
    
    [samePS(i),~,diffPS(i,2:3)] = ttest2(psH,psL);
    diffPS(i,1) = nanmean(psH)-nanmean(psL);
    if samePS(i)==0
        mrk = '.';
    else
        mrk = 'o';
    end
    plot(i,diffPS(i,1),mrk,'Color','k');
    plot([i i],diffPS(i,2:3),'k');
%     plot(i,nanmean(psL),mrk,'Color','b'); plot(i+0.2,nanmean(psH),mrk,'Color','r');
end

%% Peak Speed Flipped
figure; hold on; 
[samePSf,stdsub] = deal(zeros(length(PEAKSPEED),1));
[diffPSf] = deal(zeros(length(PEAKSPEED),3));
[PS,newLI,psF] = deal(cell(length(PEAKSPEED),1));
for i = G
    
    [okinds,stdsub(i)] = flip_dists(PEAKSPEED{G(i)});
    PS{i} = PEAKSPEED{G(i)}(okinds,:);
    
    newLI{i}{1} = okinds(PEAKSPEED{G(i)}(okinds,2)==max(PEAKSPEED{G(i)}(:,2)));
    newLI{i}{2} = okinds(PEAKSPEED{G(i)}(okinds,2)==min(PEAKSPEED{G(i)}(:,2)));
    
    lindsF = @(x) find(x(:,2)==max(x(:,2)));
    hindsF = @(x) find(x(:,2)==min(x(:,2)));
    
    linds = lindsF(PS{i});
    hinds = hindsF(PS{i});
    
    psF{i,1} = PS{i}(linds,1);
    psF{i,2} = PS{i}(hinds,1);
    
    [samePSf(i),~,diffPSf(i,2:3)] = ttest2(psF{i,2},psF{i,1});
    diffPSf(i,1) = nanmean(psF{i,2})-nanmean(psF{i,1});
    if samePSf(i)==0
        mrk = '.';
    else
        mrk = 'o';
    end
    plot(i,diffPSf(i,1),mrk,'Color','k');
    plot([i i],diffPSf(i,2:3),'k');
%     plot(i,nanmean(psL),mrk,'Color','b'); plot(i+0.2,nanmean(psH),mrk,'Color','r');
end

%% Peak Speed EQUAL
figure; hold on; 
[samePS,stdsub] = deal(zeros(length(PEAKSPEED),1));
diffPS = zeros(length(PEAKSPEED),2);
[PS,newLI,PSo] = deal(cell(length(PEAKSPEED),1));
for i = 1:length(PEAKSPEED)
    
    [~,~,okinds,stdsub(i)] = match_dists(PEAKSPEED{i});
    PS{i} = PEAKSPEED{i}(okinds,:);
    
    newLI{i}{1} = okinds(PEAKSPEED{i}(okinds,2)==max(PEAKSPEED{i}(:,2)));
    newLI{i}{2} = okinds(PEAKSPEED{i}(okinds,2)==min(PEAKSPEED{i}(:,2)));
    
    lindsF = @(x) find(x(:,2)==max(x(:,2)));
    hindsF = @(x) find(x(:,2)==min(x(:,2)));
    
    linds = lindsF(PS{i});
    hinds = hindsF(PS{i});
    
    PSo{i,1} = PS{i}(linds,1);
    PSo{i,2} = PS{i}(hinds,1);
    
    [samePS(i),~,diffPS(i,:)] = ttest2(PSo{i,2},PSo{i,1});
    if samePS(i)==0
        mrk = '.';
    else
        mrk = 'o';
    end
    plot(i,nanmean(PSo{i,2})-nanmean(PSo{i,1}),mrk,'Color','k');
    plot([i i],diffPS(i,:),'k');
%     plot(i,nanmean(psL),mrk,'Color','b'); plot(i+0.2,nanmean(psH),mrk,'Color','r');
end
%% Peak Speed Flipped
figure; hold on; 
[samePSf,stdsub] = deal(zeros(length(PEAKSPEED),1));
[diffPSf] = deal(zeros(length(PEAKSPEED),3));
[PS,newLI,psF] = deal(cell(length(PEAKSPEED),1));
for i = G
    
    [okinds,stdsub(i)] = flip_dists(PEAKSPEED{G(i)});
    PS{i} = PEAKSPEED{G(i)}(okinds,:);
    
    newLI{i}{1} = okinds(PEAKSPEED{G(i)}(okinds,2)==max(PEAKSPEED{G(i)}(:,2)));
    newLI{i}{2} = okinds(PEAKSPEED{G(i)}(okinds,2)==min(PEAKSPEED{G(i)}(:,2)));
    
    lindsF = @(x) find(x(:,2)==max(x(:,2)));
    hindsF = @(x) find(x(:,2)==min(x(:,2)));
    
    linds = lindsF(PS{i});
    hinds = hindsF(PS{i});
    
    psF{i,1} = PS{i}(linds,1);
    psF{i,2} = PS{i}(hinds,1);
    
    [samePSf(i),~,diffPSf(i,2:3)] = ttest2(psF{i,2},psF{i,1});
    diffPSf(i,1) = nanmean(psF{i,2})-nanmean(psF{i,1});
    if samePSf(i)==0
        mrk = '.';
    else
        mrk = 'o';
    end
    plot(i,diffPSf(i,1),mrk,'Color','k');
    plot([i i],diffPSf(i,2:3),'k');
%     plot(i,nanmean(psL),mrk,'Color','b'); plot(i+0.2,nanmean(psH),mrk,'Color','r');
end