kinstart = 6;
kinstop = 7;
KEEP = 1200;

speed_L = nan(length(Linds),ceil(max(1000*(comp_tt(:,kinstop)-comp_tt(:,kinstart)))));
speed_H = nan(length(Hinds),ceil(max(1000*(comp_tt(:,kinstop)-comp_tt(:,kinstart)))));

for i = 1:length(Linds)
    
    BDF_start = find(BDF.pos(:,1)>comp_tt(Linds(i),kinstart),1,'first');
    BDF_end = find(BDF.pos(:,1)>comp_tt(Linds(i),kinstop),1,'first');
    
    velx = BDF.vel(BDF_start:BDF_end,2);
    vely = BDF.vel(BDF_start:BDF_end,3);
    
    speed_L(i,1:length(velx)) = sqrt(velx.^2+vely.^2);
    
end

for i = 1:length(Hinds)
    
    BDF_start = find(BDF.pos(:,1)>comp_tt(Hinds(i),kinstart),1,'first');
    BDF_end = find(BDF.pos(:,1)>comp_tt(Hinds(i),kinstop),1,'first');
    
    velx = BDF.vel(BDF_start:BDF_end,2);
    vely = BDF.vel(BDF_start:BDF_end,3);
    
    speed_H(i,1:length(velx)) = sqrt(velx.^2+vely.^2);
    
end

speed_L = speed_L(:,1:KEEP);
speed_H = speed_H(:,1:KEEP);

av_speed_L = nanmean(speed_L,1);
av_speed_H = nanmean(speed_H,1);

std_speed_L = nanstd(speed_L,0,1);
std_speed_H = nanstd(speed_H,0,1);

figure;  hold on; 
plot(av_speed_L,'b','LineWidth',3);
plot(repmat(av_speed_L',1,2) + [std_speed_L' -std_speed_L'],'b','LineWidth',0.25);

plot(av_speed_H,'r','LineWidth',3);
plot(repmat(av_speed_H',1,2) + [std_speed_H' -std_speed_H'],'r','LineWidth',0.25);

xlabel('Time from GO (ms)','FontSize',14);
ylabel('Speed (cm/s)','FontSize',14);
title('Mr. T (3/25/2013)','FontSize',16);

%%    

for i = 1:length(comp_tt)
    
    bdf_end_ind = find(BDF.pos(:,1)<comp_tt(i,kinstop),1,'last');
    
    end_pos = BDF.pos(bdf_end_ind,2:3); % final position
    comp_tt(i,10) = atan2(end_pos(2),end_pos(1));
    
end

%% Peak Speed per Direction/Trial

max_speed_L = max(speed_L,[],2);
max_speed_H = max(speed_H,[],2);

max_all(Linds) = max_speed_L;
max_all(Hinds) = max_speed_H;

figure; plot(comp_tt(Linds,10),max_speed_L,'b.');
hold on; plot(comp_tt(Hinds,10),max_speed_H,'r.');
xlabel('Estimate Location'); ylabel('Peak Speed');

figure; plot(Linds,max_speed_L,'b.');
hold on; plot(Hinds,max_speed_H,'r.');
xlabel('Trial'); ylabel('Peak Speed');
