tt = alldays(2).tt;

%% reach direction, time to target, and max speed
%[speeds,xp,yp,xv,yv] = kin_exam(alldays(1).bdfM,tt,[],[],[6 7]);

reach_dirs = tt(:,10);
ttt = tt(:,7)-tt(:,6);
%max_speed = cellfun(@max,speeds);

KINS = [reach_dirs, ttt];%, max_speed];

%% Gain
FD = FR_change{1};

%% Width
PPC; close;
% near_inds = find(abs(circ_dist(sthet*pi/180,0))<=pi/4);
% far_inds = find(abs(circ_dist(sthet*pi/180,0))>=3*pi/4);
% 
% av_near = mean(X(:,near_inds),2);
%% av_far = mean(X(:,far_inds),2);
WE = nan(size(X,1),1);
[val_pk,allps] = deal(zeros(size(X,1),2));
for q = 1:size(X,1)
    XINPUT = X(q,:);
    recruitment_interp;
    
    allps(q,:) = p;
    
    val_pk(q,:) = [.5+p(1)-0.5*p(2) , .5+p(1)+0.5*p(2)];
    
    WE(q) = p(1)./p(2);
end

%WE = av_far - av_near;


%% Timing
Av_spike_time;

TI = nanmean(timing_array,2);

%%
tags = tt(:,3);
tagset = unique(tags);

XS = { KINS , 'kinematics'; 
      [KINS FD], 'kin + FD';
      [KINS WE], 'kin + Width';
      [KINS TI], 'kin + Timing';
      [KINS FD WE], 'kin + FD + Width';
      [KINS FD TI], 'kin + FD + Timing';
      [KINS WE TI], 'kin + Width + Timing';
      [KINS FD WE TI], 'all';
      %[KINS dfl],'dfl';
      [tt(:,11)],'visual'};
     
[bounds,preds,posts,cortag] = deal(cell(size(XS,1),1));
AUC = zeros(size(XS,1),length(tagset));
figure; hold on; 
for i = 1:size(XS,1)

    [preds{i},posts{i},cortag{i}] = naive_Bayes_cv(XS{i,1},tags,[]); % Do leave one out

    [bounds{i}(1), bounds{i}(2)] = boot_bounds(1000,@mean,preds{i}==tags,2.5,97.5);

    targs = zeros(length(tagset),length(tags));
    for j = 1:length(tagset)
        targs(j,tags==tagset(j)) = 1;
    end

%     [tpr,fpr] = roc(targs,posts{i}');
%     for j = 1:length(tpr)
%         tp = tpr{j};
%         fp = fpr{j};
%         
%         AUC(i,j) = sum(diff(fp).*0.5.*(tp(1:(end-1))+tp(2:end)));
%     end
    
    plot([i,i],[bounds{i}],'LineWidth',2);
    
    plot(i,mean(preds{i}==tags),'.','MarkerSize',5);
    
    xlim([0.5 size(XS,1)+0.5]);
end
set(gca,'XTick',1:size(XS,1),'XTickLabel',XS(:,2));



