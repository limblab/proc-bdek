kin = cds.kin;

kin.x = kin.x-2;
kin.y = kin.y+34.5;
typelist = cds.trials.cursShift+0.1*cds.trials.tgtShift;
trialtypes = unique(typelist);

typeinds = cell(length(trialtypes),1);
for i = 1:length(trialtypes)
    typeinds{i} = find(typelist==trialtypes(i) & ismember(cds.trials.result,'R') & ~isnan(cds.trials.goTime));
end
%%
dP2L = @(P1,P2,Z) (abs((P2(2)-P1(2))*Z(:,1)-(P2(1)-P1(1))*Z(:,2)+P2(1)*P1(2) - ...
                  P2(2)*P1(1))./sqrt((P2(2)-P1(2))^2 + (P2(1)-P1(1))^2)).*...
                  sign(sum(cross([repmat(P2-P1,size(Z,1),1) zeros(size(Z,1),1)],...
                                 [Z-repmat(P1,size(Z,1),1)  zeros(size(Z,1),1)],2),2));
[targloc,ctloc,deviation,TR,maxvel] = deal(cell(length(trialtypes),1));
shift_ang = 0;
offvecs = cell(length(trialtypes),1);
for i = 1:length(trialtypes)
    
    curtt = cds.trials{typeinds{i},[2 3 7 8 9 10 11]};
    
    for j = 1:size(curtt,1)
    
        starttime = find(kin.t < curtt(j,3),1,'last');
        fintime = find(kin.t < curtt(j,2),1,'last');
        
        trace_h = [kin.x(starttime:fintime) kin.y(starttime:fintime)];
        trace_c = trace_h - repmat(2*curtt(j,6)*[cos(shift_ang) sin(shift_ang)],size(trace_h,1),1);
        radius = sqrt(sum((trace_h-repmat(trace_h(1,:),size(trace_h,1),1)).^2,2));
        speed = sqrt(sum(kin.vx(starttime:fintime).^2 + kin.vy(starttime:fintime).^2,2));
        onset = find(diff(speed) > 0.25,1,'first');
        offset = find(radius > 0.8*max(radius),1,'first');
        
        TR{i}{j} = trace_h;%(radius < 6,:);

        ctloc{i}(j,:) = [0 0] + 2*curtt(j,7)*[cos(shift_ang) sin(shift_ang)];
        targloc{i}(j,:) = 6*[cos(deg2rad(curtt(j,5))) sin(deg2rad(curtt(j,5)))] + 2*curtt(j,7)*[cos(shift_ang) sin(shift_ang)];
        
%         offvec = dP2L(ctloc{i}(j,:),targloc{i}(j,:),trace_c);
%         offvec = dP2L(trace_c(onset,:),targloc{i}(j,:),trace_c(onset:offset,:));
        offvec = dP2L(trace_c(onset,:),trace_c(offset,:),trace_c(onset:offset,:));
        offvecs{i}{j,1} = offvec;
        offvecs{i}{j,2} = curtt(j,5);
        if onset > offset
            deviation{i}(j,:) = [nan nan];
            maxvel{i}(j,:) = nan;
        else
            deviation{i}(j,:) = [offvec(abs(offvec) == max(abs(offvec))) curtt(j,5)];
            maxvel{i}(j,:) = [max(speed) curtt(j,5)];
        end
    end
end

%%
nullcond = typeinds{3};
nulltt = cds.trials{nullcond,[2 3 7 8 9 10 11]};
nulldev = deviation{3}(:,1);
nullspd = maxvel{3}(:,1);

dirs = unique(nulltt(:,5));
devbydir = zeros(length(dirs),2);
spdbydir = zeros(length(dirs),2);
for i = 1:length(dirs)
    devbydir(i,:) = [dirs(i) nanmean(nulldev(nulltt(:,5)==dirs(i)))];
    spdbydir(i,:) = [dirs(i) nanmean(nullspd(nulltt(:,5)==dirs(i)))];
end
[adjdev,adjspd] = deal(cell(length(deviation),1));
for i = 1:length(deviation)
    curtt = cds.trials{typeinds{i},[2 3 7 8 9 10 11]};
    
    for j = 1:size(curtt,1);
        tloc = find(devbydir(:,1)==curtt(j,5));
        adjdev{i}(j,:) = [curtt(j,5) deviation{i}(j)-devbydir(tloc,2)];
        adjspd{i}(j,:) = [curtt(j,5) maxvel{i}(j)-spdbydir(tloc,2)];
    end
end
        
%%
[dir_dev,dir_spd] = deal(cell(length(adjdev),1));
for i = 1:length(adjdev) % trial type
    
    for j = 1:length(dirs) % directions
        
        dir_dev{i}(j,:) = [dirs(j), nanmean(adjdev{i}(adjdev{i}(:,1)==dirs(j),2))];
        dir_spd{i}(j,:) = [dirs(j), nanmean(adjspd{i}(adjspd{i}(:,1)==dirs(j),2))];
    end
end

%%
TRnrm = cell(length(TR),2);
AVtrc = cell(length(TR),1);
tdr = unique(cds.trials.tgtDir);
for i = 1:length(TR)
    for j = 1:length(TR{i})
        Xp = interp1(1:size(TR{i}{j},1),TR{i}{j}(:,1),linspace(1,size(TR{i}{j},1),100));
        Yp = interp1(1:size(TR{i}{j},1),TR{i}{j}(:,2),linspace(1,size(TR{i}{j},1),100));
        
        TRnrm{i,1}(j,:) = Xp - Xp(1);
        TRnrm{i,2}(j,:) = Yp - Yp(1);
    end
end
for i = 1:size(TR,1)
    cur_drs = cds.trials{typeinds{i},9};
    for j = 1:length(tdr)
        dinds = find(cur_drs == tdr(j));
        AVtrc{i}{j} = [nanmean(TRnrm{i,1}(dinds,:))', nanmean(TRnrm{i,2}(dinds,:))'];

    end
end
