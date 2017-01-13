%%
prior_dir = 90;

speed = nan(length(ttE),ceil(max(1000*(ttE(:,7)-ttE(:,6)))));
posx = nan(length(ttE),ceil(max(1000*(ttE(:,7)-ttE(:,6)))));
posy = nan(length(ttE),ceil(max(1000*(ttE(:,7)-ttE(:,6)))));

for i = 1:length(ttE)
    
    bdfstart = find(BDFE.bdf.pos(:,1) > ttE(i,6),1,'first');
    bdfend = find(BDFE.bdf.pos(:,1) > ttE(i,7),1,'first');
    
    velx = BDFE.bdf.vel(bdfstart:bdfend,2);
    vely = BDFE.bdf.vel(bdfstart:bdfend,3);
    
    speed_prof = sqrt(velx.^2 + vely.^2)';
    posx(i,1:length(speed_prof)) = BDFE.bdf.pos(bdfstart:bdfend,2)';
    posy(i,1:length(speed_prof)) = BDFE.bdf.pos(bdfstart:bdfend,3)';
    
    speed(i,1:length(speed_prof)) = speed_prof;
    
    maxind = find(speed(i,1:500) == max(speed(i,1:500)));
    minind = find(speed(i,500:750) == min(speed(i,500:750)));
    
    maxE.x1(i) = posx(i,maxind);
    maxE.y1(i) = posy(i,maxind);
    maxE.x2(i) = posx(i,maxind) + velx(maxind);
    maxE.y2(i) = posy(i,maxind) + vely(maxind);
    maxE.velx(i) = velx(maxind);
    maxE.vely(i) = vely(maxind);
    maxE.min(i) = speed(i,minind+499);
    
    %quiver(maxE.x(i),maxE.y(i),maxE.velx(i),maxE.vely(i));
    
end

cols = {'b' 'r' 'g' 'k'};
[IDX, C, sumd, D]  = kmeans([maxE.x1' maxE.y1' maxE.velx' maxE.vely'],3);
clust_dirs = atan2(C(:,3),C(:,2))*180./pi;
best_clust = find(abs(clust_dirs - prior_dir) == min(abs(clust_dirs - prior_dir)));
keepinds = find(IDX == best_clust);

best_ds = D(keepinds,best_clust);
CLUSTER.mean = mean(best_ds);
CLUSTER.lim = max(best_ds);
CLUSTER.centroid = C(best_clust,:);

figure; hold on;
for i = 1:length(ttE)
     plot([maxE.x1(i) maxE.x2(i)],[maxE.y1(i) maxE.y2(i)],'-o','Color',cols{IDX(i)});
end

figure; hold on;
for i = 1:length(keepinds)
    quiver(maxE.x1(keepinds(i)),maxE.y1(keepinds(i)),maxE.velx(keepinds(i)),maxE.vely(keepinds(i)));
end

avspeed = nanmean(speed(keepinds,:),1);
stdspeed = nanstd(speed(keepinds,:),1);

ttE_keep = ttE(keepinds,:);


%%
pspeed = nan(length(ttP),ceil(max(1000*(ttP(:,7)-ttP(:,6)))));
pposx = nan(length(ttP),ceil(max(1000*(ttP(:,7)-ttP(:,6)))));
pposy = nan(length(ttP),ceil(max(1000*(ttP(:,7)-ttP(:,6)))));
dist_from_cent = zeros(length(ttP),1);
for i = 1:length(ttP)
    
    bdfstart = find(BDFP.bdf.pos(:,1) > ttP(i,6),1,'first');
    bdfend = find(BDFP.bdf.pos(:,1) > ttP(i,7),1,'first');
    
    velx = BDFP.bdf.vel(bdfstart:bdfend,2);
    vely = BDFP.bdf.vel(bdfstart:bdfend,3);
    
    speed_prof = sqrt(velx.^2 + vely.^2)';
    pposx(i,1:length(speed_prof)) = BDFP.bdf.pos(bdfstart:bdfend,2)';
    pposy(i,1:length(speed_prof)) = BDFP.bdf.pos(bdfstart:bdfend,3)';
    
    pspeed(i,1:length(speed_prof)) = speed_prof;
    
    maxind = find(pspeed(i,1:500) == max(pspeed(i,1:500)));
    minind = find(pspeed(i,500:750) == min(pspeed(i,500:750)));
    
    maxP.x1(i) = pposx(i,maxind);
    maxP.y1(i) = pposy(i,maxind);
    maxP.x2(i) = pposx(i,maxind) + velx(maxind);
    maxP.y2(i) = pposy(i,maxind) + vely(maxind);
    maxP.velx(i) = velx(maxind);
    maxP.vely(i) = vely(maxind);
    maxP.min(i) = pspeed(i,minind+499);
    
    dist_from_cent(i) = sum((CLUSTER.centroid - ...
        [maxP.x1(i) maxP.y1(i) maxP.velx(i) maxP.vely(i)]).^2);
    
end

keep_prior_trials = find(dist_from_cent < CLUSTER.lim);

ttP_keep = ttP(keep_prior_trials,:);

figure; hold on;
for i = 1:length(ttP)
     plot([maxP.x1(i) maxP.x2(i)],[maxP.y1(i) maxP.y2(i)],'-o','Color',cols{ismember(i,keep_prior_trials)+1});
end

figure; hold on;
for i = 1:length(keepinds)
    quiver(maxP.x1(keepinds(i)),maxP.y1(keepinds(i)),maxP.velx(keepinds(i)),maxP.vely(keepinds(i)));
end
    
    
pavspeed = nanmean(pspeed(keep_prior_trials,:),1);
pstdspeed = nanstd(pspeed(keep_prior_trials,:),1);

figure; hold on;
plot(avspeed,'LineWidth',2); 
plot(avspeed + stdspeed,'LineWidth',0.5);
plot(avspeed - stdspeed,'LineWidth',0.5);

plot(pavspeed,'g','LineWidth',2); 
plot(pavspeed + pstdspeed,'g','LineWidth',0.5);
plot(pavspeed - pstdspeed,'g','LineWidth',0.5);

xlim([0 750]);

    
    

