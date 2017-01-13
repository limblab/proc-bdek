FIRING = horzcat(firing_absolute{1}{:});%horzcat(Firings{:});
FIRING_CO = horzcat(firing_absolute{1}{:});

tt_fire = alldays(1).tt;

unit_cell = struct2cell(vertcat(neurons{end}{:}));    
tuning_array = vertcat(unit_cell{1,:});
% 
modulations = (max(tuning_array,[],2) - min(tuning_array,[],2))';

cofire = vertcat(FIRING_CO{:});

goods = find(modulations > MODthresh & mean(cofire)>1);

EI = horzcat(EorI{:});
GOOD = cell(length(FIRING),1);
for bin = 1:length(FIRING)
    
      good_EI = find(EI(:,bin)>=0);
      GOOD{bin} = goods(ismember(goods,good_EI));
      %GOOD{bin} = E;
end

%% Likelihoods
likes = flipud(unique(tt_fire(:,3)));
for i = 1:length(likes)
    like_ind{1}{i} = find(tt_fire(:,3)==likes(i));
end
%% Variables
[fire_all,dirs_all,Fire_dir,FIRE_dir] = deal(cell(length(likes),1));

spatial_cents = -pi:pi/32:pi;
space_size = diff(spatial_cents(1:2));
[VMparams,VMparamsL,VMparamsH,curveboot,fire_trial,trlCurve,trlBase,trlGain,trlWidth,...
    outL,outH,outT,BM,GM,WM,B,G,W,C] = ...
    deal(cell(length(likes),1));
%% Run
for lik = 1:length(likes)
    
    inds = find(tt_fire(:,3)==likes(lik));
    
    for bin = 1:length(FIRING)

        dPD = dfrompd(:,GOOD{bin});%;{1}{end};
        dirs_all{lik}{bin} = reshape(dPD(inds,:),[],1);

        FIRE_norm_good = (FIRING{bin}(inds,GOOD{bin}))./...
                repmat(var(cofire(:,GOOD{bin})),length(inds),1);


        fire_all{lik}{bin}  = reshape(FIRE_norm_good,[],1);

        fire_trial{lik}{bin} = FIRE_norm_good;

        clc; fprintf('fitting...\ncondition: %d/%d\ntime bin: %d/%d\n',lik,length(likes),bin,length(FIRING));

    end
    
    for dir = 1:length(spatial_cents)
        for bin = 1:length(FIRING)
            dirinds = find(abs(circ_dist(dirs_all{lik}{bin},spatial_cents(dir)))<(space_size/2));
            Fire_dir{lik}(dir,bin) = nanmean(fire_all{lik}{bin}(dirinds),1);
            %FIRE_dir{lik}{dir} = fire_all{lik}(dirinds,:);
        end
    end 
end

%%
c2p = {'k','b','r','r'};
mincol = min(min(horzcat(Fire_dir{:})));
maxcol = max(max(horzcat(Fire_dir{:})));

time_cell = cellfun(@(x) (x(1:end-1)+x(2:end))/2,tune_ranges,'UniformOutput',0);
time_cell{2} = time_cell{2}+1000;
time_cents = horzcat(time_cell{:});

figure; hold on; 
for i = 1:length(likes)
    
    subplot(1,length(likes),i); 
    imagesc(time_cents,spatial_cents,Fire_dir{i},[mincol,maxcol]);
    
    hold on;
    plot([0 0],[-1 1]*(pi+pi/64),'w','LineWidth',3);
    patch([800 1000 1000 800],[-1 -1 1 1]*(pi+pi/32),'w','EdgeColor','w');
    box off
end
% %% Diff from baseline activity
% % find baseline
% baseinds = find(time_cents<0);
% basemeans = nanmean(cell2mat(cellfun(@(x) x(:,baseinds),Fire_dir,'UniformOutput',0)'),2);
% baserep = repmat(basemeans,1,size(Fire_dir{1},2));
% 
% Fire_dir_diff = cell(size(Fire_dir));
% for i = 1:length(likes)
%     Fire_dir_diff{i} = Fire_dir{i}-baserep;
% end
% 
% mincol2 = min(min(horzcat(Fire_dir_diff{:})));
% maxcol2 = max(max(horzcat(Fire_dir_diff{:})));
% 
% figure; hold on; 
% for i = 1:length(likes)
%     
%     subplot(1,length(likes),i); 
%     imagesc(time_cents,spatial_cents,Fire_dir_diff{i},[mincol2,maxcol2]);
%     
%     hold on;
%     plot([0 0],[-1 1]*(pi+pi/64),'w','LineWidth',3);
%     patch([800 1000 1000 800],[-1 -1 1 1]*(pi+pi/32),'w','EdgeColor','w');
%     box off
% end
%%

if 0

figure; hold on;
facea = .25;
edgea = 0;
% BASELINES = cellfun(@(x) x(:,1),VMparams,'UniformOutput',0); BASELINES = horzcat(BASELINES{:});
% GAINS = cellfun(@(x) x(:,2),VMparams,'UniformOutput',0); GAINS = horzcat(GAINS{:});
% WIDTHS = cellfun(@(x) x(:,3),VMparams,'UniformOutput',0); WIDTHS = horzcat(WIDTHS{:});
subplot(1,3,1); hold on; title('Baseline','FontSize',18);
for i = 1:length(likes)
    plot(time_cents,BM{i},c2p{i}); 
    patch([time_cents fliplr(time_cents)],[B{i}(1,:) fliplr(B{i}(2,:))],c2p{i},'FaceAlpha',facea,'EdgeAlpha',edgea);
    yl = ylim;
    plot([0 0],[yl(1) yl(2)],'w','LineWidth',3);
    patch([800 1000 1000 800],[yl(1) yl(1) yl(2) yl(2)],'w','EdgeColor','w');
    xlim([time_cents(1) time_cents(end)]);
end

subplot(1,3,2); hold on; title('Gain','FontSize',18);
for i = 1:length(likes) 
    plot(time_cents,GM{i},c2p{i}); 
    patch([time_cents fliplr(time_cents)],[G{i}(1,:) fliplr(G{i}(2,:))],c2p{i},'FaceAlpha',facea,'EdgeAlpha',edgea);
    yl = ylim;
    plot([0 0],[yl(1) yl(2)],'w','LineWidth',3);
    patch([800 1000 1000 800],[yl(1) yl(1) yl(2) yl(2)],'w','EdgeColor','w');
    xlim([time_cents(1) time_cents(end)]);
end

subplot(1,3,3); hold on; title('Width','FontSize',18);
for i = 1:length(likes) 
    plot(time_cents,WM{i},c2p{i}); 
    patch([time_cents fliplr(time_cents)],[W{i}(1,:) fliplr(W{i}(2,:))],c2p{i},'FaceAlpha',facea,'EdgeAlpha',edgea);
    yl = ylim;
    plot([0 0],[yl(1) yl(2)],'w','LineWidth',3);
    patch([800 1000 1000 800],[yl(1) yl(1) yl(2) yl(2)],'w','FaceAlpha',1,'EdgeColor','w');
    xlim([time_cents(1) time_cents(end)]);
end
    
%%

lastind = find(time_cell{1}<800,1,'last');
% indbegin = 7;
indbegin = 3;

[bs,ppcfit,XCOR,firereshape,bs_cos] = deal(cell(length(likes),1));
[R2_VM,R2_cos] = deal(zeros(length(likes),1));
figure; hold on; 
for i = 1:length(likes)
    
%     xsfire = repmat(spatial_cents',8,1);
%     ysfire = reshape(Fire_dir{i}(:,(lastind-7):lastind),[],1);
    
    xsfire = spatial_cents';
    ysfire = mean(Fire_dir{i}(:,(lastind-indbegin):lastind),2);
    
    firecell = cellfun(@(x) x(:,(lastind-indbegin):lastind),FIRE_dir{i},'UniformOutput',0);
    
    firereshape{i} = cellfun(@(x) reshape(x,[],1),firecell,'UniformOutput',0);
 
    xrep = cell(length(firereshape),1);
    for j = 1:length(firereshape{i})
        xrep{j} = ones(size(firereshape{i}{j})).*spatial_cents(j);
    end
    
    Firecat = vertcat(firereshape{i}{:});
    xcat = vertcat(xrep{:});
    
    %%% Cosine (exp) Fitting
%     helperfunc = @(x) cos_helper_func(x(:,1),x(:,2));
%     
%     [xl,xu,xboot] = boot_bounds(100, helperfunc, [xcat,Firecat],2.5,97.5);
%    
%     firetime = FIRE_dir{i}((lastind-7):lastind);
%     YSFire = vertcat(firetime{:});
 
    %%% VM Fitting
    vm_func = @(xs,p) p(1) + p(2)*exp(p(3)*cos(xs-p(4))); 
%     bs_cos{i} = glmfit([cos(xsfire) sin(xsfire)],ysfire,'poisson');
%     bs{i} = VM_fit(xsfire,ysfire);
    bs{i} = VM_fit(xcat,Firecat);
    
    ppcfit{i} = vm_func(-pi:0.01:pi,bs{i});
    
    %[xl,xu,xboot] = boot_bounds(100, @VM_helper_func, [xcat,Firecat],2.5,97.5);
%     [psl,psu] = boot_bounds(100, @VM_helper_func_kap, [xcat, Firecat],2.5,97.5);
    
    %ppcfit{i} = glmval(bs{i},[cos(-pi:0.01:pi)' sin(-pi:0.01:pi)'],'log');
    
    plot(xsfire,ysfire,'.','Color',c2p{i});
    plot(-pi:0.01:pi,ppcfit{i},c2p{i},'LineWidth',3);
    
    SSE = sum((ysfire-vm_func(xsfire,bs{i})).^2);
    SST = sum((ysfire-mean(ysfire)).^2);
    
%     SSE_cos = sum((ysfire-glmval(bs_cos{i},[cos(xsfire) sin(xsfire)],'log')).^2);
    
    R2_VM(i) = 1 - SSE/SST;
%     R2_cos(i) = 1 - SSE_cos/SST;
    
    %patch([-pi:0.01:pi pi:-0.01:-pi],[xl' fliplr(xu')],c2p{i});
    
%     for j = 1:length(spatial_cents)
%        %xshift = repmat(circshift(spatial_cents',j),8,1);
%        yshift = reshape(circshift(Fire_dir{i}(:,(lastind-7):lastind),j),[],1);
%        
%        XCOR{i}(j) = circ_corrcc(yshift,ysfire);
%        
%     end
%     XCOR{i} = [1 XCOR{i}];
%     
%     plot(XCOR{i},c2p{i});
end
   
end
