% pdi = 2;
TT = alldays(pdi).tt;
if isfield(alldays,'kin'); KIN = alldays(1).kin; else KIN = alldays(1).bdfM; end
% cop = {'b','g','r'};
% likes = flipud(unique(TT(:,3)));

% figure; hold on;
% speed = cell(length(likes),1);
% V = zeros(3,length(likes));
% for li = 1:length(likes)
% 
%     l_inds = find(TT(:,3)==likes(li)); 
%     
%     for tri = 1:length(l_inds)
%         
%         i1 = find(KIN.pos(:,1) > TT(l_inds(tri),6),1,'first');
%         i2 = find(KIN.pos(:,1) < TT(l_inds(tri),7),1,'last');
%         
%         xvel = KIN.vel(i1:i2,2)';
%         yvel = KIN.vel(i1:i2,3)';
%         
%         speed{li}(tri) = max(sqrt(xvel.^2 + yvel.^2));
%         clc; fprintf('%d/%d - %d/%d\n',li,length(likes),tri,length(l_inds));
%     end
% %     n = histc(speed{li},0:0.5:50);
% %     bar(0:0.5:50,n,cop{li});
% 
%     V(1,li) = nanmean(speed{li});
%     [V(2,li), V(3,li)] = boot_bounds(1000,@nanmean,speed{li},2.5,97.5);
% 
% end
[react_time,topspeed, time_topspeed] = deal(nan(size(TT,1),1));
POS = [smooth(KIN.pos(:,2),50) smooth(KIN.pos(:,3),50)];

poscent = zeros(size(alldays(1).tt,1),2);
for tri = 1:size(alldays(1).tt,1)
    i0 = find(KIN.pos(:,1) > alldays(1).tt(tri,5),1,'first');
    i1 = find(KIN.pos(:,1) > alldays(1).tt(tri,9),1,'first');

    poscent(tri,:) = mean(POS(i0:i1,:),1);
end
POS = POS - repmat(mean(poscent,1),size(POS,1),1);

% figure; hold on; 
for tri = 1:size(TT,1)

    i0 = find(KIN.pos(:,1) > TT(tri,5),1,'first');
    i1 = find(KIN.pos(:,1) > TT(tri,9),1,'first');
    i2 = find(KIN.pos(:,1) < TT(tri,10),1,'last');

    xvel = smooth(KIN.vel(i1:i2,2),50)';
    yvel = smooth(KIN.vel(i1:i2,3),50)';
    
    xpos = POS(i1:i2,1)';
    ypos = POS(i1:i2,2)';
    
    xrest = POS(i0:i1,1)';
    yrest = POS(i0:i1,2)';
    
    xall = POS(i0:i2,1)';
    yall = POS(i0:i2,2)';
    
    radis = sqrt(xpos.^2+ypos.^2);
    radis_all = sqrt(POS(i0:i2,1).^2+POS(i0:i2,2).^2);
%     radis2 = sqrt((xpos-xpos(1)).^2+(ypos-ypos(1)).^2);
%     radis0 = sqrt((xrest-xrest(1)).^2+(yrest-yrest(1)).^2);
%     radisall = sqrt(xall.^2+yall.^2);
%     
%     dradisall = diff(radisall);
%     dradisrest = dradisall(1:(i1-i0-1));
%     
%     ddradisall = smooth(diff(dradisall),50);
%     ddradisrest = ddradisall(50:500);
%     
%     ddradthresh = 3*max(abs(ddradisrest));
    
    spd = sqrt(xvel.^2 + yvel.^2);
    
%     rest_radis(tri) = max(sqrt(xrest.^2+yrest.^2));
    
%     pos_as = atan2(ypos,xpos);
%     pos2_as = atan2(ypos-ypos(1),xpos-xpos(1));
%     vel_as = atan2(yvel,xvel);
    
%     errang = circ_dist(pos_as,alldays(2).tt(tri,10));
%     errend = sqrt((xpos-xpos(end)).^2 + (ypos-ypos(end)).^2);
    
%     acc_time(tri) = find(ddradisall>ddradthresh,1,'first')-(i1-i0);

    react_time(tri) = find(radis_all>1.4,1,'first')-(i1-i0);
    topspeed(tri) = max(spd);
    time_topspeed(tri) = find(spd == max(spd),1,'first');
%     clc; fprintf('%d/%d - %d/%d\n',li,length(likes),tri,length(l_inds));
%     plot(radis); pause; cla;
    clc; fprintf('%d/%d\n',tri,length(TT));
end