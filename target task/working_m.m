typeCOL{1} = 'b';
typeCOL{2} = 'r';
typeCOL{3} = 'k';

reaches = comp_tt(:,7)-comp_tt(:,6);
maxreach = max(ceil(1000*(reaches)));

waits = comp_tt(:,6)-comp_tt(:,5);
maxwait = max(ceil(1000*(waits)));

speeds = zeros(length(comp_tt),maxreach);
waitspeeds = zeros(length(comp_tt),maxwait);

vxs = zeros(length(comp_tt),maxreach);
vys = zeros(length(comp_tt),maxreach);

vxafter = zeros(length(comp_tt),100);

types = zeros(length(comp_tt),1);

targpos = zeros(length(comp_tt),2);
startpos = zeros(length(comp_tt),2);
figure; hold on;
for i = 1:length(OUTPUTS.y)
    
    BDFtarg = comp_tt(i,5);
    BDFstart = comp_tt(i,6);
    BDFend = comp_tt(i,7);
    
    bti = find(BDF.pos(:,1)>BDFtarg,1,'first');
    bsi = find(BDF.pos(:,1)>BDFstart,1,'first');
    bei = find(BDF.pos(:,1)<BDFend,1,'last');
    
    speed = sqrt(BDF.vel(bsi:bei,2).^2 + BDF.vel(bsi:bei,3).^2);
    waitspeed = sqrt(BDF.vel(bti:bsi,2).^2 + BDF.vel(bti:bsi,3).^2);
    vx = BDF.vel(bsi:bei,2);
    vy = BDF.vel(bsi:bei,3);
    vxafter(i,:) = BDF.vel(bsi:bsi+99,2)';
    
    vxs(i,1:length(vx)) = vx;
    vys(i,1:length(vy)) = vy;
    
    targpos(i,:) = BDF.pos(bti,2:3);
    startpos(i,:) = BDF.pos(bsi,2:3);
    
    speeds(i,1:length(speed)) = speed;
    waitspeeds(i,1:length(waitspeed)) = waitspeed;
    
    if ismember(i,nospikeinds)
        type = 3;
    elseif OUTPUTS.x(i,end)==1
        type = 1;
    elseif OUTPUTS.x(i,end)==2
        type = 2;
    end
    
    types(i)= type;
    
    plot(mean(vxafter(i,:)),type,'.','Color',typeCOL{type},'MarkerSize',4);
    
end

meanspeedL = mean(speeds(types==1,:),1);
    timeL = reaches(types==1);
meanspeedH = mean(speeds(types==2,:),1);
    timeH = reaches(types==2);
meanspeedN = mean(speeds(types==3,:),1);
    timeN = reaches(types==3);

meanwspeedL = mean(waitspeeds(types==1,:),1);
meanwspeedH = mean(waitspeeds(types==2,:),1);
meanwspeedN = mean(waitspeeds(types==3,:),1);

figure; hold on;
plot(meanspeedL,'b');
plot(meanspeedH,'r');
plot(meanspeedN,'k');
title('Movement');

figure; hold on;
plot(meanwspeedL,'b');
plot(meanwspeedH,'r');
plot(meanwspeedN,'k');
title('Wait');

figure; hold on;
plot(targpos(types==1,1),targpos(types==1,2),'b.','MarkerSize',3);
plot(targpos(types==2,1),targpos(types==2,2),'r.','MarkerSize',3);
plot(targpos(types==3,1),targpos(types==3,2),'k.','MarkerSize',3);
title('Position at Target');

figure; hold on;
plot(startpos(types==1,1),startpos(types==1,2),'b.','MarkerSize',3);
plot(startpos(types==2,1),startpos(types==2,2),'r.','MarkerSize',3);
plot(startpos(types==3,1),startpos(types==3,2),'k.','MarkerSize',3);
title('Position at Go');

meanvxL = mean(vxs(types==1,:),1);
meanvxH = mean(vxs(types==2,:),1);
meanvxN = mean(vxs(types==3,:),1);

meanvyL = mean(vys(types==1,:),1);
meanvyH = mean(vys(types==2,:),1);
meanvyN = mean(vys(types==3,:),1);

figure; 
subplot(1,2,1); hold on;

plot(meanvxL,'b');
plot(meanvxH,'r');
plot(meanvxN,'k');

subplot(1,2,2); hold on;
plot(meanvyL,'b');
plot(meanvyH,'r');
plot(meanvyN,'k');

title('Vel')
