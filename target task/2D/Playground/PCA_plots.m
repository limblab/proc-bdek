%% Find 'stay' and 'swap' trials
area2use = 'PMd_prep';

pretarg = 1:8;
init = 8:49;
precue = 8:63;
delay1 = 8:28;
delay2 = 29:49;
delay = 8:49;
mem = 49:63;
move = 103:111;
postcue = 64:size(M.(lower(area2use)){2}.d,2);
timeall = 1:size(M.(lower(area2use)){2}.d,2);
timetarg = 8:size(M.(lower(area2use)){2}.d,2);

region = 30:38;%mem;
trials = [];

align_col = 10;
[targids] = deal(cell(length(M.tt),1));
for i = 1:length(M.tt)
    targids{i} = mod(round(mod(M.tt{i}(:,align_col)+4*pi,2*pi)./(pi/4))+1,9);
    targids{i}(targids{i}==0)= 1;
end

for i = 1:length(M.(lower(area2use)))
 
    trials.stays{i}.all = find(reshape(nanmean(diff(M.(lower(area2use)){i}.d([1 5],region,:),1),2),[],1)>0);
    trials.swaps{i}.all = find(reshape(nanmean(diff(M.(lower(area2use)){i}.d([1 5],region,:),1),2),[],1)<0);

%     trials.corrects{i}.all = find(diff(floor(M.tt{i}(:,18:19)),[],2) == 0);
%     trials.wrongs{i}.all = find(diff(floor(M.tt{i}(:,18:19)),[],2) ~= 0);
%     trials.all{i}.all = [trials.corrects{i}.all; trials.wrongs{i}.all];
    
     trials.corrects{i}.all = find(diff(floor(M.tt{i}(:,9:10)),[],2) == 0);
    trials.wrongs{i}.all = find(diff(floor(M.tt{i}(:,9:10)),[],2) ~= 0);
    trials.all{i}.all = [trials.corrects{i}.all; trials.wrongs{i}.all];
    for j = 1:8

        trials.targs{i}.targs{j,:} = find(targids{i}==j);

        trials.stays{i}.targs{j,:} = trials.stays{i}.all(targids{i}(trials.stays{i}.all)==j);
        trials.swaps{i}.targs{j,:} = trials.swaps{i}.all(targids{i}(trials.swaps{i}.all)==j);

        trials.corrects{i}.targs{j,:} = trials.corrects{i}.all(targids{i}(trials.corrects{i}.all)==j);
        trials.wrongs{i}.targs{j,:} = trials.wrongs{i}.all(targids{i}(trials.wrongs{i}.all)==j);
        trials.all{i}.targs{j,:} = [trials.corrects{i}.targs{j,:}; trials.wrongs{i}.targs{j,:}];
    end
    
    for j = 1:8
        trials.prefer{i}.targs{j,:} = [trials.stays{i}.targs{j}; -trials.swaps{i}.targs{mod(j+3,8)+1}];
    end
end
CLENS = M.CLENS;
LENS = M.LENS;
%% Plot targ A and B representations
block = 2 ;
areas = {'pmd_prep','pmd_move'};%,'m1'};
trials2use = L;%'corrects';%biases.notany{2}.all;%'corrects';%'corrects';%reps.BnotA{block};%'swaps';%notanybias;%'swaps';%'all';%'swaps';%strongbiasA(3);%'corrects';%'all';%;'corrects';%'corrects';%'swaps'; %'prefer';% wrongs{4}];%stayswap{2,6};
target2use = 'all';
timeper = timeall;
roundto = .1;
reference_targs = [3 7];%[2 3 4 6 7];

[lb, ub] = deal(cell(2,1));
figure; hold on;
cap = @(x,limit) min(abs(x),limit).*sign(x);
round2nearest = @(x,n) n*ceil(x./n);
for i = 1:length(areas)

    subplot(1,length(areas),i); hold on; 
    title(areas{i},'FontSize',18);
    if strcmp(target2use,'all'); target2use = [1 2 3 4 5 6 7 8]; end
    
    if isnumeric(trials2use)
        pertarg = cellfun(@(x) x(ismember(x,abs(trials2use))),trials.all{block}.targs,'Uni',0);
        alltrialsign = cellfun(@(x) sign(trials2use(ismember(abs(trials2use),x(ismember(x,abs(trials2use)))))),trials.all{block}.targs,'Uni',0);
        trialsign = cell2mat(alltrialsign(abs(target2use)));
        trialset = abs(cell2mat(pertarg(abs(target2use))));
        targsign = cell2mat(cellfun(@(x) repmat(x(1),x(2),1),mat2cell([sign(target2use); ...
                            cellfun(@(x) size(x,1),pertarg(abs(target2use)))'],2,ones(length(target2use),1)),...
                            'UniformOutput',0)').*trialsign;

    else
        trialset = abs(cell2mat(trials.(trials2use){block}.targs(abs(target2use))));
        signlist = sign(cell2mat(trials.(trials2use){block}.targs(abs(target2use))));

        targsign = cell2mat(cellfun(@(x) repmat(x(1),x(2),1),mat2cell([sign(target2use); ...
                            cellfun(@(x) size(x,1),trials.(trials2use){block}.targs(abs(target2use)))'],2,ones(length(target2use),1)),...
                            'UniformOutput',0)').*signlist;
    end

    blockmat = M.(areas{i}){block}.d(:,timeper,trialset);
    blockmat(:,:,targsign<0) = circshift(blockmat(:,:,targsign<0),-4,1);
    
    T = [max(timeper(1)-5,0) min(timeper(end)+5,CLENS(end))];
    xs = timeper;

    ys = nanmean(blockmat([5 1],:,:)./repmat(nanmean(blockmat(reference_targs,:,:),1),2,1,1),3)';
%     ys = nanmean(blockmat([5 1],:,:),3)';

    plot(repmat(xs',1,2),ys,'-');
    if length(trialset)>1

%         [lb{1},ub{1}] = boot_bounds(1000,@(x) nanmean(x,1),...
%                     (squeeze(blockmat(5,:,:)) - squeeze(nanmean(blockmat(reference_targs,:,:),1)))',2.5,97.5);
%   
%         [lb{2},ub{2}] = boot_bounds(1000,@(x) nanmean(x,1),...
%                     (squeeze(blockmat(1,:,:)) - squeeze(nanmean(blockmat(reference_targs,:,:),1)))',2.5,97.5);
                
        [lb{1},ub{1}] = boot_bounds(1000,@(x) nanmean(x,1),...
                    (squeeze(blockmat(5,:,:))./squeeze(nanmean(blockmat(reference_targs,:,:),1)))',2.5,97.5);
  
        [lb{2},ub{2}] = boot_bounds(1000,@(x) nanmean(x,1),...
                    (squeeze(blockmat(1,:,:))./squeeze(nanmean(blockmat(reference_targs,:,:),1)))',2.5,97.5);
                
%         [lb{1},ub{1}] = boot_bounds(1000,@(x) nanmean(x,1),...
%                     (squeeze(blockmat(5,:,:)))',2.5,97.5);
%   
%         [lb{2},ub{2}] = boot_bounds(1000,@(x) nanmean(x,1),...
%                     (squeeze(blockmat(1,:,:)))',2.5,97.5);

        Ylimit = round2nearest(max(reshape(vertcat([lb{:}; ub{:}]),[],1)),roundto);
        Ylimitl = -round2nearest(max(reshape(vertcat(-[lb{:}; ub{:}]),[],1)),roundto);
        
        plot([xs fliplr(xs)],cap([lb{1}' fliplr(ub{1}')],Ylimit),'k');
        plot([xs fliplr(xs)],cap([lb{2}' fliplr(ub{2}')],Ylimit),'k');

    elseif isfield(M.(areas{i}){block},'AB')
        lbub1 = squeeze(M.(areas{i}){block}.AB(1,:,trialset,:));
        lbub2 = squeeze(M.(areas{i}){block}.AB(2,:,trialset,:));
        Ylimit = round2nearest(max(reshape([lbub1;lbub2],[],1)),roundto);
        Ylimitl = -round2nearest(max(reshape(-[lbub1;lbub2],[],1)),roundto);
        plot([xs fliplr(xs)],cap([lbub1(:,1)' fliplr(lbub1(:,2)')],Ylimit),'k');
        plot([xs fliplr(xs)],cap([lbub2(:,1)' fliplr(lbub2(:,2)')],Ylimit),'k');
    else
        Ylimit = round2nearest(max(reshape(ys,[],1)),roundto);
        Ylimitl = -round2nearest(max(reshape(-ys,[],1)),roundto);
    end
    

%     plot(repmat(cumsum(LENS),2,1),repmat([-Ylimit;Ylimit],1,length(LENS)),'Color',[.5 .5 .5]);%
    plot(repmat(cumsum(LENS),2,1),repmat([Ylimitl;Ylimit],1,length(LENS)),'Color',[.5 .5 .5]);
    
    plot([T(1) T(2)],[1 1],'k-','LineWidth',0.1)
    xlim([T(1) T(2)])
%     ylim([-Ylimit Ylimit]);
    ylim([Ylimitl Ylimit]);
end
%% Plot diff between targ A and B 
block = 2;
areas = {'pmd5','m15'};%,'pmdbad'};%,'m1r','m1c'};
trials2use = 'corrects';% wrongs{4}];%stayswap{2,6};
target2use = {'all'};%{[1 5],[2 6],[3 7],[4 8]};%{1,2,3,4,5,6,7,8};%{1,2,3,4,5,6,7,8};
show_bounds = 'on';
timeper = timeall;

% trialsets = {corrects{2}(1),wrongs{2}(2)};
ys = cell(length(target2use),length(areas));
% [lb, ub] = deal(cell(length(trialsets),1));
figure; hold on;
for i = 1:length(areas)

    subplot(1,length(areas),i); hold on; 
    title(areas{i},'FontSize',18); 
    
    T = [max(timeper(1)-5,0) min(timeper(end)+5,CLENS(end))];
    xs = timeper;

    for j = 1:length(target2use)
        if strcmp(target2use{j},'all')
            target2use{j} = 1:8;
        end
        
        if isnumeric(trials2use)
            trialset = abs(trials2use);
            trialsign = sign(trials2use);
        else
            trialset = cell2mat(trials.(trials2use){block}.targs(abs(target2use{j})));

            trialsign = cell2mat(cellfun(@(x) repmat(x(1),x(2),1),mat2cell([sign(target2use{j}); ...
                                cellfun(@(x) size(x,1),trials.(trials2use){block}.targs(abs(target2use{j})))'],2,ones(length(target2use{j}),1)),...
                                'UniformOutput',0)');
        end
                        
        blockmat = M.(areas{i}){block}.d(:,timeper,trialset);
        blockmat(:,:,trialsign<0) = circshift(blockmat(:,:,trialsign<0),-4,1);
        
        if size(blockmat,3) > 1
            ys{j,i} = nanmean(squeeze(diff(blockmat([1 5],:,:),[],1)),2);
        else
            ys{j,i} = diff(blockmat([1 5],:,:),[],1);
        end
        plot(xs,ys{j,i},'-');
        if length(trialset)>1 && ~strcmp(show_bounds,'off')
            [lb,ub] = boot_bounds(100,@(x) nanmean(x,1),squeeze(diff(blockmat([1 5],:,:),[],1))',2.5,97.5);
            plot([xs fliplr(xs)],[lb' fliplr(ub')],'k');
        end
    end
    yl = getfield(gca,'YLim');%1.1*abs(max(ys(:)));
    plot(repmat(cumsum(LENS),2,1),repmat([yl(1);yl(2)],1,length(LENS)),'Color',[.5 .5 .5])
    plot([T(1) T(2)],[0 0],'k-','LineWidth',0.1)
    xlim([T(1) T(2)])
    ylim(yl)
end

%% Plot targ A and B representations for Planning v Execution
block = 2;
areas = {'pmd','m1'};
reptypes = {1,2};
trials2use = 'corrects';% wrongs{4}];%stayswap{2,6};
target2use = 4;

pltstyl = {'-','.-'};
[lb, ub] = deal(cell(2,1));
ys = cell(length(reptypes),length(areas));
subplotorder = reshape(reshape(1:(length(areas)*length(reptypes)),length(areas),length(reptypes))',1,[]);
figure; hold on;
for i = 1:length(areas)
    for j = 1:length(reptypes)
        subplot(length(reptypes),length(areas),subplotorder(length(reptypes)*(i-1)+j)); hold on; 
        title(sprintf('%s - epoch:%d',areas{i},reptypes{j}),'FontSize',18);

        blockmat = M.(areas{i}){block}.('ep'){reptypes{j}}; 
        
        xs = 1:size(blockmat,2);
        trialset = trials.(trials2use){block}.targs{target2use};
        ys{j,i} = nanmean(blockmat([5 1],:,trialset),3)' - repmat(nanmean(nanmean(blockmat([3 7],:,trialset),1),3)',1,2);
        plot(ys{j,i},pltstyl{j});
        
        if length(trialset)>1
            [lb{1},ub{1}] = boot_bounds(100,@(x) nanmean(x,1),...
                        (reshape(blockmat(5,:,trialset),size(blockmat,2),length(trialset))...
                      - reshape(nanmean(blockmat([3 7],:,trialset),1),size(blockmat,2),length(trialset)))',2.5,97.5);
            [lb{2},ub{2}] = boot_bounds(100,@(x) nanmean(x,1),...
                        (reshape(blockmat(1,:,trialset),size(blockmat,2),length(trialset))...
                      - reshape(nanmean(blockmat([3 7],:,trialset),1),size(blockmat,2),length(trialset)))',2.5,97.5);

            plot([xs fliplr(xs)],[lb{1}' fliplr(ub{1}')],'k');
            plot([xs fliplr(xs)],[lb{2}' fliplr(ub{2}')],'k');
        end
        
        plot(repmat(cumsum(LENS),2,1),repmat([-0.05;.15],1,length(LENS)),'Color',[.5 .5 .5])
        plot([0 max(CLENS)],[0 0],'k-','LineWidth',0.1)
    end
% xlim([0 126])
% ylim([-0.05 0.15])
end

%% Plot both PMd and M1
yd = cellfun(@(x) -diff(x,[],2),ys,'UniformOutput',0);
yds{1} = [yd{1,1}, yd{2,1}]; yds{2} = [yd{1,2} yd{2,2}];
absmax = @(A,dim) max(abs(A),[],dim).*(2*(max(abs(A),[],dim)==max(A,[],dim))-1);

AvB = cellfun(@(x) absmax(x,2),yds,'UniformOutput',0);
PE = cellfun(@(x) sum(abs(x),2)',ys,'UniformOutput',0);
for i = 1:size(PE,2); PEc{i} = vertcat(PE{:,i})'; end
% PvE = cellfun(@(x) (x(:,2)-x(:,1))./sum(x,2),PEc,'UniformOutput',0);
PvE = cellfun(@(x) x(:,2)-x(:,1),PEc,'UniformOutput',0);

%
xP = AvB{1};
yP = PvE{1};
xM = AvB{2};
yM = PvE{2};


figure; hold on; subplot(1,2,1); hold on; title('PMd');
xb1 = 1.1*max(abs(xP));
yb1 = 1.1*max(abs(yP));
xlim([-1 1]*xb1);
ylim([-1 1]*yb1); 
hAx1 = gca; set(hAx1,'XTickLabel','','YTickLabel','',...
                     'Xtick',[],'Ytick',[],'xcolor','w','ycolor','w');

subplot(1,2,2); hold on; title('M1');
xb2 = 1.1*max(abs(xM));
yb2 = 1.1*max(abs(yM));
xlim([-1 1]*xb2); 
ylim([-1 1]*yb2); 
hAx2 = gca; set(hAx2,'XTickLabel','','YTickLabel','',...
                     'Xtick',[],'Ytick',[],'xcolor','w','ycolor','w');
plot(hAx1,[0 0],[-1 1]*yb1,'k'); plot(hAx1,[-1 1]*xb1,[0 0],'k');
% text(hAx1,.1,-yb1*1.1,'\leftarrow planning','FontSize',16);
% text(hAx1,0.55,-yb1*1.1,'execution \rightarrow','FontSize',16)
% ha = text(hAx1,-0.1,yb1/2,'A \rightarrow','FontSize',16); set(ha,'rotation',90);
% hb = text(hAx1,-0.1,-yb1/2,'\leftarrow B','FontSize',16); set(hb,'rotation',90);
plot(hAx2,[0 0],[-1 1]*yb2,'k'); plot(hAx2,[-1 1]*xb2,[0 0],'k');
% text(hAx2,.1,-yb2*1.1,'\leftarrow planning','FontSize',16);
% text(hAx2,0.55,-yb2*1.1,'execution \rightarrow','FontSize',16)
for i = 2:126
    if i > 111
        clr = 'r';
    elseif i > 103
        clr = 'm';
    elseif i > 63
        clr = 'g';
    elseif i > 49
        clr = 'k';
    elseif i > 8
        clr = 'b';
    else
        clr = 'y';
    end

    plot(hAx1,xP((i-1):i),yP((i-1):i),'.-','Color',clr,'MarkerSize',20); 
%     plot(hAx1,nanmean(PE_state(tr,i)),nanmean(plotmet(tr,i)),'.','Color',clr,'MarkerSize',nanmean(PE_weight(1,i))*40);
    plot(hAx2,xM((i-1):i),yM((i-1):i),'.-','Color',clr,'MarkerSize',20); 

%     plot(hAx2,nanmean(PE_stateM(tr,i)),nanmean(plotmetM(tr,i)),'.','Color',clr,'MarkerSize',nanmean(PE_weightM(1,i))*40);

    pause(0.1); 
end    

%% Plot both PMd and M1 both targets
block = 2;
areas = {'pmd','m1'};
trialset = wrongs{2};%wrongs{2}(2);%corrects{2};%stayswap{2,6};

ys = cell(length(reptypes),length(areas));
for i = 1:length(areas)
    for j = 1:2
        areamats = M.(areas{i});
        blockmat = M.(areas{i}){block}.('ep'){j}; 
        ys{j,i} = nanmean(blockmat([5 1],:,trialset),3)' - repmat(nanmean(nanmean(blockmat([3 7],:,trialset),1),3)',1,2);
    end
end
%
xP = ys{1,1};
yP = ys{2,1};
xM = ys{1,2};
yM = ys{2,2};

figure; hold on; subplot(1,2,1); hold on; title('PMd');
xb1 = 1.1*max(abs(xP(:)));
yb1 = 1.1*max(abs(yP(:)));
xlim([-1 1]*xb1);
ylim([-1 1]*yb1); 
hAx1 = gca; set(hAx1,'XTickLabel','','YTickLabel','',...
                     'Xtick',[],'Ytick',[],'xcolor','w','ycolor','w');

subplot(1,2,2); hold on; title('M1');
xb2 = 1.1*max(abs(xM(:)));
yb2 = 1.1*max(abs(yM(:)));
xlim([-1 1]*xb2); 
ylim([-1 1]*yb2); 
hAx2 = gca; set(hAx2,'XTickLabel','','YTickLabel','',...
                     'Xtick',[],'Ytick',[],'xcolor','w','ycolor','w');
plot(hAx1,[0 0],[-1 1]*yb1,'k'); plot(hAx1,[-1 1]*xb1,[0 0],'k');
% text(hAx1,.1,-yb1*1.1,'\leftarrow planning','FontSize',16);
% text(hAx1,0.55,-yb1*1.1,'execution \rightarrow','FontSize',16)
% ha = text(hAx1,-0.1,yb1/2,'A \rightarrow','FontSize',16); set(ha,'rotation',90);
% hb = text(hAx1,-0.1,-yb1/2,'\leftarrow B','FontSize',16); set(hb,'rotation',90);
plot(hAx2,[0 0],[-1 1]*yb2,'k'); plot(hAx2,[-1 1]*xb2,[0 0],'k');
% text(hAx2,.1,-yb2*1.1,'\leftarrow planning','FontSize',16);
% text(hAx2,0.55,-yb2*1.1,'execution \rightarrow','FontSize',16)
clrs = {'b','k','r','g','g','y','y'};
c = [0, 0.447, 0.741; 0.85, 0.325, 0.098];
for i = 2:120

    if ismember(i,CLENS)
        plot(hAx1,xP(i,1),yP(i,1),'.','Color',clrs{i==CLENS},'MarkerSize',50);
        plot(hAx1,xP(i,2),yP(i,2),'.','Color',clrs{i==CLENS},'MarkerSize',50);

        plot(hAx2,xM(i,1),yM(i,1),'.','Color',clrs{i==CLENS},'MarkerSize',50); 
        plot(hAx2,xM(i,2),yM(i,2),'.','Color',clrs{i==CLENS},'MarkerSize',50); 
    end
        
    plot(hAx1,xP((i-1):i,1),yP((i-1):i,1),'.-','Color',c(1,:),'MarkerSize',20);
    plot(hAx1,xP((i-1):i,2),yP((i-1):i,2),'.-','Color',c(2,:),'MarkerSize',20);

    plot(hAx2,xM((i-1):i,1),yM((i-1):i,1),'.-','Color',c(1,:),'MarkerSize',20); 
    plot(hAx2,xM((i-1):i,2),yM((i-1):i,2),'.-','Color',c(2,:),'MarkerSize',20); 

    pause(0.1); 
end    

%% Plot both PMd and M1 A v B
block = 2;
areas = {'pmd','m1'};
trialset = wrongs{2};%wrongs{2}(2);%corrects{2};%stayswap{2,6};

ys = cell(length(reptypes),length(areas));
for i = 1:length(areas)
    for j = 1:2
        areamats = M.(areas{i});
        blockmat = M.(areas{i}){block}.('ep'){j}; 
        ys{j,i} = nanmean(blockmat([5 1],:,trialset),3)' - repmat(nanmean(nanmean(blockmat([3 7],:,trialset),1),3)',1,2);
    end
end
%
xP = -absmax([diff(ys{1,1},[],2), diff(ys{2,1},[],2)],2);
xM = -absmax([diff(ys{1,2},[],2), diff(ys{2,2},[],2)],2);
yP = absmax(ys{2,1},2) - absmax(ys{1,1},2);
yM = absmax(ys{2,2},2) - absmax(ys{1,2},2);

figure; hold on; subplot(1,2,1); hold on; title('PMd');
xb1 = 1.1*max(abs(xP(:)));
yb1 = 1.1*max(abs(yP(:)));
xlim([-1 1]*xb1);
ylim([-1 1]*yb1); 
hAx1 = gca; set(hAx1,'XTickLabel','','YTickLabel','',...
                     'Xtick',[],'Ytick',[],'xcolor','w','ycolor','w');

subplot(1,2,2); hold on; title('M1');
xb2 = 1.1*max(abs(xM(:)));
yb2 = 1.1*max(abs(yM(:)));
xlim([-1 1]*xb2); 
ylim([-1 1]*yb2); 
hAx2 = gca; set(hAx2,'XTickLabel','','YTickLabel','',...
                     'Xtick',[],'Ytick',[],'xcolor','w','ycolor','w');
plot(hAx1,[0 0],[-1 1]*yb1,'k'); plot(hAx1,[-1 1]*xb1,[0 0],'k');
% text(hAx1,.1,-yb1*1.1,'\leftarrow planning','FontSize',16);
% text(hAx1,0.55,-yb1*1.1,'execution \rightarrow','FontSize',16)
% ha = text(hAx1,-0.1,yb1/2,'A \rightarrow','FontSize',16); set(ha,'rotation',90);
% hb = text(hAx1,-0.1,-yb1/2,'\leftarrow B','FontSize',16); set(hb,'rotation',90);
plot(hAx2,[0 0],[-1 1]*yb2,'k'); plot(hAx2,[-1 1]*xb2,[0 0],'k');
% text(hAx2,.1,-yb2*1.1,'\leftarrow planning','FontSize',16);
% text(hAx2,0.55,-yb2*1.1,'execution \rightarrow','FontSize',16)
clrs = {'b','k','r','m','g','y','y'};
c = [0, 0.447, 0.741; 0.85, 0.325, 0.098];
for i = 2:CLENS(end-1)

    if ismember(i,CLENS)
        plot(hAx1,xP(i,1),yP(i,1),'.','Color',clrs{i==CLENS},'MarkerSize',50);

        plot(hAx2,xM(i,1),yM(i,1),'.','Color',clrs{i==CLENS},'MarkerSize',50); 
    end
        
    plot(hAx1,xP((i-1):i,1),yP((i-1):i,1),'.-','Color',c(1,:),'MarkerSize',10);

    plot(hAx2,xM((i-1):i,1),yM((i-1):i,1),'.-','Color',c(1,:),'MarkerSize',10); 

    pause(0.1); 
end    


%% Plot timepoint estimates, PMd and M1
blocks = {1,2};
areas = {'pmd','m1'};%,'m1'};
trials2use = 'corrects';% wrongs{4}];%stayswap{2,6};
target2use = 'all';

ys = cell(length(blocks),length(areas));
figure; hold on;
cap = @(x,limit) min(abs(x),limit).*sign(x);
for i = 1:length(areas)

    subplot(1,length(areas),i); hold on; 
    title(areas{i},'FontSize',18);
    
    for bl = 1:length(blocks)
        blockmat = M.(areas{i}){blocks{bl}}.t;

        if strcmp(target2use,'all'); target2use = 1:8; end
        trialset = cell2mat(trials.(trials2use){blocks{bl}}.targs(target2use));
        T = size(blockmat,2);
        xs = 1:T;
        ys{bl,i} = nanmean(blockmat(trialset,:));
        plot(repmat(xs',1,2),ys{bl,i},'-');
        if length(trialset)>1
            [lb,ub] = boot_bounds(100,@(x) nanmean(x,1),blockmat(trialset,:),2.5,97.5);

            Ylimit = 10*ceil(max(abs(reshape(vertcat([lb; ub]),[],1)))./10);

            plot([xs fliplr(xs)],cap([lb' fliplr(ub')],Ylimit),'k');
        end
    end
    plot([0 max(CLENS)],[0 max(CLENS)],'k');
    plot(repmat(cumsum(LENS),2,1),repmat([0;T],1,length(LENS)),'Color',[.5 .5 .5])

    plot([0 T],[0 0],'k-','LineWidth',0.1)
    xlim([0 T])
    ylim([0 T]);
%     ylim([-0.05 0.15])
% ylim([-0.2 0.25])
end
%%
figure; hold on;  
for i = 1:126; 
    cla; 
    plot(ys{1,1}(i),-.25,'b.','MarkerSize',50); 
    plot(ys{2,1}(i),.25,'r.','MarkerSize',50); 
    ylim([-1 1]); 
    plot(repmat(CLENS(2:end),2,1),repmat([-1; 1],1,length(LENS)),'k'); 
    plot([0 126],[0 0],'k','LineWidth',0.25); 
    pause(0.05); 
end


%% Plot diff between targ A and B multtargs
block = 2;
areas = {'pmd','m1'};
trials2use = 'corrects';% wrongs{4}];%stayswap{2,6};
target2use = {[1 5],[2 6],[3 7],[4 8]};
targflip = [1, -1];

% trialsets = {corrects{2}(1),wrongs{2}(2)};
flipinds = cell(length(target2use{1}),1);
ys = cell(length(target2use),length(areas));
[lb, ub] = deal(cell(length(trialsets),1));
figure; hold on;
for i = 1:length(areas)
    subplot(1,length(areas),i); hold on; 
    title(areas{i},'FontSize',18); 

    xs = 1:size(blockmat,2);

    for j = 1:length(target2use)
        trialinds = trials.(trials2use){block}.targs(target2use{j});
        for fl = 1:length(targflip)
            flipinds{fl} = repmat(targflip(fl),length(trialinds{fl}),size(M.(areas{i}){block}.d,2));
        end
        blockmat = reshape(diff(M.(areas{i}){block}.d([1 5],:,cell2mat(trialinds)),[],1),[],size(cell2mat(trialinds),1))'.*cell2mat(flipinds);
 
        ys{j,i} = nanmean(blockmat)';
        plot(xs,ys{j,i},'-');
        if length(trialset)>1
            [lb,ub] = boot_bounds(100,@(x) nanmean(x,1),blockmat,2.5,97.5);
            plot([xs fliplr(xs)],[lb' fliplr(ub')],'k');
        end
    end
    yl = getfield(gca,'YLim');%1.1*abs(max(ys(:)));
    plot(repmat(cumsum(LENS),2,1),repmat([yl(1);yl(2)],1,length(LENS)),'Color',[.5 .5 .5])
    plot([0 max(xs)],[0 0],'k-','LineWidth',0.1)
    xlim([0 max(xs)])
    ylim(yl)
end



%% Plot diff between targ A and B 
block = 1;
areas = {'pmd'};%,'m1'};%,'m1'};%,'pmdbad'};%,'m1r','m1c'};
trialblocks2use = {shortlull,longlull,aftercatch};%'corrects'};%{biasAboth,biasBboth};%{biasAatcue_m1, biasBatcue_m1,nobiasatcue_m1};%{trials.wrongs{2}.targs{3}};%{biases.anyA{block}.all,biases.anyB{block}.all};%{biases.strongA{block}.all,biases.weakA{block}.all,biases.strongB{block}.all,biases.weakB{block}.all,biases.notany{block}.all};%{'corrects'};%{biases.anyA{block},biases.notany{block},biases.anyB{block}};%, 
% trialblocks2use = {};
target2use = {'all'};%[1],[2],[3],[4],[5],[6],[7],[8]};%,[6]};%,[2 6],[3 7],[4 8]};%{1,2,3,4,5,6,7,8};%{1,2,3,4,5,6,7,8};
show_bounds = 'off';
timeper = timeall;

% trialsets = {corrects{2}(1),wrongs{2}(2)};
ys = cell(length(target2use),length(areas));
% [lb, ub] = deal(cell(length(trialsets),1));
figure; hold on;
for i = 1:length(areas)

    subplot(1,length(areas),i); hold on; 
    title(areas{i},'FontSize',18); 
    
    T = [max(timeper(1)-5,0) min(timeper(end)+5,CLENS(end))];
    xs = timeper;
    
    for k = 1:length(trialblocks2use)
        
        trials2use = trialblocks2use{k};

        for j = 1:length(target2use)
            if strcmp(target2use{j},'all')
                target2use{j} = 1:8;
            end

            if isnumeric(trials2use)
                pertarg = cellfun(@(x) x(ismember(x,abs(trials2use))),trials.all{block}.targs,'Uni',0);
                alltrialsign = cellfun(@(x) sign(trials2use(ismember(abs(trials2use),x(ismember(x,abs(trials2use)))))),trials.all{block}.targs,'Uni',0);
                trialsign = cell2mat(alltrialsign(abs(target2use{j})));
                trialset = abs(cell2mat(pertarg(abs(target2use{j}))));
                targsign = cell2mat(cellfun(@(x) repmat(x(1),x(2),1),mat2cell([sign(target2use{j}); ...
                            cellfun(@(x) size(x,1),pertarg(abs(target2use{j})))'],2,ones(length(target2use{j}),1)),...
                            'UniformOutput',0)').*trialsign;

%                 trialset = abs(trials2use);
%                 trialsign = sign(trials2use);
            else
                trialset = cell2mat(trials.(trials2use){block}.targs(abs(target2use{j})));

                targsign = cell2mat(cellfun(@(x) repmat(x(1),x(2),1),mat2cell([sign(target2use{j}); ...
                                    cellfun(@(x) size(x,1),trials.(trials2use){block}.targs(abs(target2use{j})))'],2,ones(length(target2use{j}),1)),...
                                    'UniformOutput',0)');
            end

            blockmat = M.(areas{i}){block}.d(:,timeper,trialset);
            blockmat(:,:,targsign<0) = circshift(blockmat(:,:,targsign<0),-4,1);
            targnums = targids{block}(trialset);
%             biascomp = zeros(length(targnums),length(timeper));
%             for tn = 1:length(targnums)
%                 biascomp(tn,:) = daxis(targnums(tn),timeper);
%             end
            
            if size(blockmat,3) > 1
                ys{j,i} = nanmean(squeeze(diff(blockmat([1 5],:,:),[],1)),2);
%                 ys{j,i} = nanmean(squeeze(diff(blockmat([1 5],:,:),[],1))-biascomp',2);
            else
                ys{j,i} = diff(blockmat([1 5],:,:),[],1);
            end
            if ~isempty(ys{j,i})
                plot(xs,ys{j,i},'-');
                if length(trialset)>1 && ~strcmp(show_bounds,'off')
                    [lb,ub] = boot_bounds(100,@(x) nanmean(x,1),squeeze(diff(blockmat([1 5],:,:),[],1))',2.5,97.5);
                    plot([xs fliplr(xs)],[lb' fliplr(ub')],'k');
                end
            else
                plot(xs,nan(size(xs)),'-');
            end
        end
    end
    yl = getfield(gca,'YLim');%1.1*abs(max(ys(:)));
    plot(repmat(cumsum(LENS),2,1),repmat([yl(1);yl(2)],1,length(LENS)),'Color',[.5 .5 .5])
    plot([T(1) T(2)],[0 0],'k-','LineWidth',0.1)
    xlim([T(1) T(2)])
    ylim(yl)
end


