%%
ttO = tt;
posi = @(p) round(polyval([1000 -999],p)); %function for turning times into indexes 

% Fill out position, velocity and acceleration variables
pos = bdf.pos; pos(:,2) = pos(:,2)-3; pos(:,3) = pos(:,3)+33;% correct the xy offsets 
vel = bdf.vel; 
acc = bdf.acc;
if isfield(bdf,'force'); forc = bdf.force; else forc = []; end

% 
for i = 1:size(tt,1)
    ps = pos(posi(tt(i,9)),2:3)-pos(posi(tt(i,8)),2:3);
    moveangs(i) = atan2(ps(2),ps(1));
    
    %tt(i,13) = atan2(ps(2),ps(1));
end
%bump.out = .601;


% Initialize variables
[bump_ang,bump_mag,bump_frc,bump_frcx,bump_frcy] = deal(zeros(size(tt,1),1));
[bumpsx,bumpsy,bvelx,bvely,bspeed,baccx,baccy,maxspeed_dir,bforcex,bforcey,forcemags] = ...
    deal(nan(size(tt,1),round(bump.out*1000)+1));
 
figure; hold on; 
for i = 1:size(tt,1)

    t1 = tt(i,7);%-0.2; % time at bump onset
    t2 = t1 + bump.out; % time at bump offset
    
    % Fill out hand position at all times during the bump
    bumpsx(i,:) = (pos(posi(t1):posi(t2),2)-pos(posi(t1),2))';
    bumpsy(i,:) = (pos(posi(t1):posi(t2),3)-pos(posi(t1),3))';
    % Plot the bump trajectory
    plot(bumpsx(i,:),bumpsy(i,:));
    
    
    % Center the trajectory on xy position at onset
    xbump = pos(posi(t2),2)-pos(posi(t1),2);
    ybump = pos(posi(t2),3)-pos(posi(t1),3);
    
    % Fill out velocity, acceleration, and speed trajectories during bump
    bvelx(i,:) = vel(posi(t1):posi(t2),2);
    bvely(i,:) = vel(posi(t1):posi(t2),3);
    
    baccx(i,:) = acc(posi(t1):posi(t2),2);
    baccy(i,:) = acc(posi(t1):posi(t2),3);
    
    if size(forc,1)>1
        bforcex(i,:) = forc(posi(t1):posi(t2),2);
        bforcey(i,:) = forc(posi(t1):posi(t2),3);
        
        plot(bforcex(i,:),bforcey(i,:),'g');
        
        forcemags(i,:) = sqrt(bforcex(i,:).^2 + bforcey(i,:).^2);
        magforceloc = find(forcemags(i,:)==max(forcemags(i,:)));

        bump_frc(i) = max(forcemags(i,:));
        bump_frcx(i) = bforcex(magforceloc);
        bump_frcy(i) = bforcey(magforceloc);
    end

    bspeed(i,:) = sqrt(bvelx(i,:).^2 + bvely(i,:).^2);
    
    % Find point of maximum speed during bump
    maxspeed_i = find(bspeed(i,:)==max(bspeed(i,:)),1,'first');
    maxspeed_dir(i) = atan2(bvely(i,maxspeed_i),bvelx(i,maxspeed_i));
    
    %bump_ang(i) = atan2(ybump,xbump);
    bump_ang(i) = maxspeed_dir(i);
    bump_mag(i) = sqrt(ybump.^2 + xbump.^2);

    % Plot intended bump directions in red
    %plot([0 5*cos(tt(i,2))],[0, 5*sin(tt(i,2))],'r');
end
%% Separate conditions and plot behavior   
%ANGLES = maxspeed_dir;
ANGLES = bump_ang;

vis = find(tt(:,3)-floor(tt(:,3))~=0); % Trials with visual target
novis = find(tt(:,3)-floor(tt(:,3))==0); % Trials without visual target

bumpsizes = unique(tt(novis,4)); % Find how many bump sizes were used
[bumps, errors_bump] = deal(cell(length(bumpsizes),1));
for i = 1:length(bumpsizes)
    bumps{i} = novis(tt(novis,4)==bumpsizes(i)); % separate bump sizes
end
allbumps = vertcat(bumps{:}); % Create a variable with all bump sizes


% Plot endpoint location against target location (or bump direction)
figure; hold on;
subplot(1,length(bumpsizes)+1,1); hold on; 
plot(tt(vis,2)/pi*180,(mod(tt(vis,13)+10*pi,2*pi))/pi*180,'r.'); title('Visual','FontSize',16);
xlabel('Target Location','FontSize',14); ylabel('Endpoint','FontSize',14);
xlim([0 360]); ylim([0 360]);

for i= 1:length(bumpsizes)
    subplot(1,length(bumpsizes)+1,i+1); hold on; 
    plot((mod(ANGLES(bumps{i})+10*pi,2*pi))/pi*180,(mod(tt(bumps{i},13)+10*pi,2*pi))/pi*180,'b.');
    title(sprintf('Bumps (%.1f cm)',bumpsizes(i)),'FontSize',16);
    xlabel('Bump Direction','FontSize',14); ylabel('Endpoint','FontSize',14);
    xlim([0 360]); ylim([0 360]);
end

%% Categorize directions
targs = kmeans(tt(:,2),length(unique(tt(:,2))),'start',unique(tt(:,2))); % cluster our targets to find indices for each

%% Histograms
allerrors = zeros(size(tt,1),1);
errors_vis = circ_dist(tt(vis,13),tt(vis,2));
allerrors(vis) = errors_vis;

% Plot the errors as a function of direction
figure; hold on; 
subplot(2,length(bumpsizes)+1,1); hold on;
plot(mod(tt(vis,2)+10*pi,2*pi),errors_vis,'.');
ylim([-5 5]);
xlim([0 2*pi]);
for i = 1:length(bumpsizes)
    errors_bump{i} = circ_dist(tt(bumps{i},13),ANGLES(bumps{i}));
    subplot(2,length(bumpsizes)+1,i+1); hold on; 
    plot(mod(ANGLES(bumps{i})+10*pi,2*pi),errors_bump{i},'b.');
    ylim([-5 5]);
    xlim([0 2*pi]);
    
    allerrors(bumps{i}) = errors_bump{i};
end

%% Calculate mean and concentration (kappa) of errors with bootstrapped confidence bounds
vis_loc = nan(length(unique(targs)),1);
bump_loc = nan(length(unique(targs)),length(bumpsizes));
[V_errs,VK_errs] = deal(nan(length(unique(targs)),2));
[B_errs,BK_errs,Eloc] = deal(cell(length(bumpsizes),1));
Vloc = nan(length(unique(targs)),3);
for i = 1:length(unique(targs))
    inds = vis(targs(vis) == i);
    
    vis_loc(i) = circ_mean(tt(inds,2));
    [V_errs(i,1),V_errs(i,2)] = boot_bounds(1000,@circ_mean,allerrors(inds),2.5,97.5);
    [VK_errs(i,1),VK_errs(i,2)] = boot_bounds(1000,@circ_kappa,allerrors(inds),2.5,97.5);
    Vloc(i,1) = circ_mean(mod(tt(inds,13)+10*pi,2*pi));
    [Vloc(i,2),Vloc(i,3)] = boot_bounds(1000,@circ_mean,mod(tt(inds,13)+10*pi,2*pi),2.5,97.5);
    
    for j = 1:length(bumpsizes)
        inds = bumps{j}(targs(bumps{j}) == i);
        
        bump_loc(i,j) = circ_mean(ANGLES(inds));
        
        [B_errs{j}(i,1),B_errs{j}(i,2)] = boot_bounds(1000,@circ_mean,allerrors(inds),2.5,97.5);
        [BK_errs{j}(i,1),BK_errs{j}(i,2)] = boot_bounds(1000,@circ_kappa,allerrors(inds),2.5,97.5);
        
        Eloc{j}(i,1) = circ_mean(tt(inds,13));
        [Eloc{j}(i,2),Eloc{j}(i,3)] = boot_bounds(1000,@circ_mean,tt(inds,13),2.5,97.5);
        
    end
end

% Plot errors as a function of direction
subplot(2,length(bumpsizes)+1,length(bumpsizes)+2);
plot(mod(repmat(vis_loc,1,2)'+10*pi,2*pi),V_errs','r','LineWidth',3);
xlim([0 2*pi]);
ylim([-1 1]);
title('Visual','FontSize',16);
xlabel('Target Direction','FontSize',14);
ylabel('Error','FontSize',14);

for i = 1:length(bumpsizes)
    
    subplot(2,length(bumpsizes)+1,i+2+length(bumpsizes));
    plot(mod(repmat(bump_loc(:,i),1,2)'+10*pi,2*pi),B_errs{i}','b','LineWidth',3);
    xlim([0 2*pi]);
    ylim([-1 1]);
    title(sprintf('Bumps (%.1f cm)',bumpsizes(i)),'FontSize',16);
    xlabel('Bump Direction','FontSize',14);
    ylabel('Error','FontSize',14);
end

%% Plot errors as a function of direction
figure; hold on;
subplot(1,length(bumpsizes)+1,1);
plot(mod(repmat(vis_loc,1,2)'+10*pi,2*pi),VK_errs','r','LineWidth',3);
xlim([0 2*pi]);
ylim([0 max(max(vertcat(BK_errs{:})))]);
title('Visual','FontSize',16);
xlabel('Target Direction','FontSize',14);
ylabel('Kappa (error)','FontSize',14);

for i = 1:length(bumpsizes)
    
    subplot(1,length(bumpsizes)+1,i+1);
    plot(mod(repmat(bump_loc(:,i),1,2)'+10*pi,2*pi),BK_errs{i}','b','LineWidth',3);
    xlim([0 2*pi]);
    ylim([0 max(max(vertcat(BK_errs{:})))]);
    title(sprintf('Bumps (%.1f cm)',bumpsizes(i)),'FontSize',16);
    xlabel('Bump Direction','FontSize',14);
    ylabel('Kappa(error)','FontSize',14);
    
end
%% Histogram of errors 
figure; hold on; 
subplot(1,length(bumpsizes)+1,1);
hist(allerrors(vis)/pi*180,50); xlim([-180 180]);
title(sprintf('Visual\nKappa: %.1f',circ_kappa(allerrors(vis))),'FontSize',16);
xlabel('Error (deg)','FontSize',14); ylabel('Count','FontSize',14);
for i = 1:length(bumpsizes)
    subplot(1,length(bumpsizes)+1,i+1);
    hist(allerrors(bumps{i})/pi*180,50); xlim([-180 180]);
    title(sprintf('Bump (%.1f cm)\nKappa: %.1f',bumpsizes(i),circ_kappa(allerrors(bumps{i}))),'FontSize',16);
    xlabel('Error (deg)','FontSize',14); ylabel('Count','FontSize',14);
end

%% 
figure; hold on; 
targs_of_interest = [2,4,6,8];
col2plot = {'r','b','g','k'};
for i = 1:length(targs_of_interest);
    toi = targs_of_interest(i);
    %nh{i} = histc(mod(tt(bumps{1}(targs(bumps{1})==toi),13)+10*pi,2*pi)/pi*180,0:10:360);
    nh{i} = histc(mod(tt(allbumps(targs(allbumps)==toi),13)+10*pi,2*pi)/pi*180,0:5:360);

    nv{i} = histc(mod(tt(vis(targs(vis)==toi),13)+10*pi,2*pi)/pi*180,0:5:360);
    
    %bar(0:10:360,nv{i},col2plot{2});
    bar(0:5:360,nh{i},col2plot{1}); hold on; 
    
    
%     plot(mod(bump_ang(bumps{1}(targs(bumps{1})==toi))+10*pi,2*pi)/pi*180,...
%          mod(tt(bumps{1}(targs(bumps{1})==toi),13)+10*pi,2*pi)/pi*180,'.','Color',col2plot{i});
end

% 
% for i = 1:length(targs_of_interest)
%     toi = targs_of_interest(i);
%     plot([1 1]*mod(circ_mean(bump_ang(bumps{1}(targs(bumps{1})==toi)))+10*pi,2*pi)/pi*180,[0 1.1*max(vertcat(nh{:}))],'--','Color',col2plot{i},'LineWidth',2);
% end
%%
figure; hold on; 
for i =1:length(bumps{1})
    
    ps = pos(posi(tt(bumps{1}(i),8)):posi(tt(bumps{1}(i),9)),2:3);
    vs = vel(posi(tt(bumps{1}(i),8)):posi(tt(bumps{1}(i),9)),2:3);
    sps = vs(:,1).^2 + vs(:,2).^2;
    msp = find(sps==max(sps));
    
    %takeoff_ang = atan2(ps(
    
    %anglemove(i) = atan2(ps(end,2)-ps(1,2),ps(end,1)-ps(1,1));
    
    plot(ps(:,1)-ps(1,1),ps(:,2)-ps(1,2),'b');
    %plot(cos(anglemove(i)),sin(anglemove(i)),'b'); 
end