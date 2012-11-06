function[bdf,tt,trials2,trains,ts_low,pva] = prepare_vars4GLM(Monkey, Day)

%PARAMS
feedback_loc = 4;
y_offset = -7;
GLOBALx = -2;
GLOBALy = 32.5;
% find the day's PMd file

if strcmp(Monkey,'MrT')
    path = '\\165.124.111.182\data\MrT_9I4\bdf';
    m = what(path); 
    
elseif strcmp(Monkey, 'Mihili')
    path = '\\165.124.111.182\data\Mihili_12A3\bdf';
    m = what(path);

else 
    fprintf('Monkey Not Found\n');
end

mats = m.mat;

findday = regexp(mats,Day);
findpmd = regexp(mats,'PMd');
findm1 = regexp(mats,'M1');

dayinds = [];
for i = 1:length(findday)
    if ~isempty(findday{i})
        dayinds = [dayinds i];
    end
end
pmdinds = [];
for i = 1:length(findpmd)
    if ~isempty(findpmd{i})
        pmdinds = [pmdinds i];
    end
end
behaveinds = [];
for i = 1:length(findm1)
    if ~isempty(findm1{i})
        behaveinds = [behaveinds i];
    end
end


PMd_ind = dayinds(ismember(dayinds,pmdinds));
bdfunit = load([path '\' mats{PMd_ind}]);

trains = spiketrains(bdfunit.bdf,1);

Behave_ind = dayinds(ismember(dayinds,behaveinds));
bdfbehav = load([path '\' mats{Behave_ind(1)}]);
bdf = bdfbehav.bdf;
% Create trial table. 
tt = getTT(bdf);
tt(isnan(tt(:,5)),:) = [];

% Create updated trial table trials2
trials2(:,1:2) = tt(:,2:3);

% Find feedback ON ts

for i = 1:length(tt)  
    start = find(bdf.pos(:,1) > tt(i,6),1,'first');
    last = find(bdf.pos(:,1) > tt(i,7),1,'first');
    mid = find(bdf.pos(start:last,3)>(-GLOBALy+y_offset+feedback_loc),1,'first')+start-1;
    
    trials2(i,3) = bdf.pos(start,1);
    trials2(i,4) = bdf.pos(mid,1);
    trials2(i,5) = bdf.pos(last,1);
end

% Create PVA
pva.pos = bdf.pos; pva.vel = bdf.vel; pva.acc = bdf.acc;

% Find minimum speed timestamps
ts_low = zeros(length(trials2),1);
for i = 1:length(tt)
    
    mid = trials2(i,4);
    last = trials2(i,5);
    
    midind = find(bdf.pos(:,1) > mid,1,'first');
    lastind = find(bdf.pos(:,1) > last,1,'first');
    
    % Speed profile
    speed = sqrt(bdf.vel(midind:lastind,2).^2 + bdf.vel(midind:lastind,3).^2);
    acc = (-speed(5:end) + 8.*speed(4:end-1) - 8.*speed(2:end-3) + speed(1:end-4))./(12*.001);
    acc = [nan; nan; acc; nan; nan];
    
    zer_acc = [];
    for j= 1:length(acc)-1
        if acc(j+1) > 0 && acc(j) < 0
            zer_acc = [zer_acc (j + midind - 1)];
            break;
        end
    end

    if length(zer_acc) == 1;
        ts_low(i) = bdf.pos(zer_acc,1);
    elseif isempty(zer_acc);
        fprintf('No Min Speed: trial %d',i);
    else
        ts_low(i) = bdf.pos(zer_accs(2),1);
        fprintf('Check Trial: %d',i);
        pause;
    end
 
end

end


    