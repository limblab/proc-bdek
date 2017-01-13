function[ydistcov] = movetraces(bdf,comp_tt)

holds = comp_tt(:,6)-comp_tt(:,5);
comp_tt(holds>1.3,:) = [];
lowinds = find(comp_tt(:,3)==min(comp_tt(:,3)));
highinds = find(comp_tt(:,3)==max(comp_tt(:,3)));

figure; hold on;

triallengths = comp_tt(:,7)-comp_tt(:,5);
longesttrial = ceil(max(triallengths)*1000)+1;
speeds = nan(length(comp_tt),longesttrial);
stds = nan(2,longesttrial);


for i = 1:length(comp_tt)
    
    subplot(1,2,1); hold on;
    tind(i) = find(bdf.pos(:,1)<comp_tt(i,5),1,'last');
    gind(i) = find(bdf.pos(:,1)<comp_tt(i,6),1,'last');
    rind(i) = find(bdf.pos(:,1)<comp_tt(i,7),1,'last');
        
    %plot(bdf.pos(tind(i):gind(i),2),bdf.pos(tind(i):gind(i),3),'b');
    plot(bdf.pos(tind(i),2),bdf.pos(tind(i),3),'bx');
    %plot(bdf.pos(gind(i):rind(i),2),bdf.pos(gind(i):rind(i),3),'g');
    plot(bdf.pos(rind(i),2),bdf.pos(rind(i),3),'kx');
    plot(bdf.pos(gind(i),2),bdf.pos(gind(i),3),'gx');
    
    speed = sqrt(bdf.vel(gind(i):rind(i),2).^2 + bdf.vel(gind(i):rind(i),3).^2);
    maxspind(i) = find(speed==max(speed));
    if ismember(i,lowinds)
        plot(bdf.pos(maxspind(i)+gind(i)-1,2),bdf.pos(maxspind(i)+gind(i)-1,3),'bo');
    elseif ismember(i,highinds)
        plot(bdf.pos(maxspind(i)+gind(i)-1,2),bdf.pos(maxspind(i)+gind(i)-1,3),'ro');
    end
    
    subplot(1,2,2); hold on;
    
    tspeed = sqrt((bdf.vel(tind(i):gind(i),2)).^2 + (bdf.vel(tind(i):gind(i),3)).^2);
    gspeed = sqrt((bdf.vel(gind(i):rind(i),2)).^2 + (bdf.vel(gind(i):rind(i),3)).^2);
    
    plot(bdf.pos(tind(i):gind(i),1)-bdf.pos(tind(i)),tspeed,'b');
    plot(bdf.pos(gind(i):rind(i),1)-bdf.pos(gind(i))+bdf.pos(gind(i),1)-bdf.pos(tind(i),1),gspeed,'g');
    
    speeds(i,1:length(tspeed)+length(gspeed)-1) = [tspeed(1:end-1);gspeed]';
    
    %pause;
    
end

Lspeed = nanmean(speeds(lowinds,:),1);
Lstd = nanstd(speeds(lowinds,:),1);

Hspeed = nanmean(speeds(highinds,:),1);
Hstd = nanstd(speeds(highinds,:),1);

figure; hold on;
plot(0:0.001:(size(Lspeed,2)-1)/1000,Lspeed,'b','LineWidth',3);
plot(0:0.001:(size(Hspeed,2)-1)/1000,Hspeed,'r','LineWidth',3);

legend('Low Uncertainty','High Uncertainty');

plot(0:0.001:(size(Lspeed,2)-1)/1000,Lspeed+Lstd,'b','LineWidth',1);
plot(0:0.001:(size(Lspeed,2)-1)/1000,Lspeed-Lstd,'b','LineWidth',1);

plot(0:0.001:(size(Hspeed,2)-1)/1000,Hspeed+Hstd,'r','LineWidth',1);
plot(0:0.001:(size(Hspeed,2)-1)/1000,Hspeed-Hstd,'r','LineWidth',1);

title('Speed Profiles','FontSize',16);
xlabel('Time from target appearance (s)','FontSize',14);
ylabel('Hand Speed (cm/s)','FontSize',14);


