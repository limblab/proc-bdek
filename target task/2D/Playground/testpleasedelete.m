function testpleasedelete(M)
%% Find 'stay' and 'swap' trials
pretarg = 1:8;
init = 8:49;
precue = 8:63;
delay1 = 8:28;
delay2 = 29:49;
delay = 8:49;
mem = 38:48;%49:63;
move = 103:111;
postcue = 64:size(M.pmd{2}.d,2);
timeall = 1:size(M.pmd{2}.d,2);
timetarg = 8:size(M.pmd{2}.d,2);

region = mem;
trials = [];
area2use = 'PMd';

align_col = 19;
[targids] = deal(cell(length(M.tt),1));
for i = 1:length(M.tt)
    targids{i} = mod(round(mod(M.tt{i}(:,align_col)+4*pi,2*pi)./(pi/4))+1,9);
    targids{i}(targids{i}==0)= 1;
end

for i = 1:length(M.(lower(area2use)))
 
    trials.stays{i}.all = find(reshape(nanmean(diff(M.(lower(area2use)){i}.d([1 5],region,:),1),2),[],1)>0);
    trials.swaps{i}.all = find(reshape(nanmean(diff(M.(lower(area2use)){i}.d([1 5],region,:),1),2),[],1)<0);

    trials.corrects{i}.all = find(diff(floor(M.tt{i}(:,18:19)),[],2) == 0);
    trials.wrongs{i}.all = find(diff(floor(M.tt{i}(:,18:19)),[],2) ~= 0);
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
LENS = M.LENS;
CLENS = M.CLENS;
%% Plot targ A and B representations
block = 2;
areas = {'pmd','m1'};
trials2use = 'all';%'corrects';%biases.notany{2}.all;%'corrects';%'corrects';%reps.BnotA{block};%'swaps';%notanybias;%'swaps';%'all';%'swaps';%strongbiasA(3);%'corrects';%'all';%;'corrects';%'corrects';%'swaps'; %'prefer';% wrongs{4}];%stayswap{2,6};
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