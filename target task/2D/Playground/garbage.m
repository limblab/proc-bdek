%% Plot diff between targ A and B 
block = 2;
areas = {'pmd','m1'};%,'m1'};%,'pmdbad'};%,'m1r','m1c'};
trialblocks2use = {trials.corrects{2}.all};%{trials.wrongs{2}.targs{3}};%{biases.anyA{block}.all,biases.anyB{block}.all};%{biases.strongA{block}.all,biases.weakA{block}.all,biases.strongB{block}.all,biases.weakB{block}.all,biases.notany{block}.all};%{'corrects'};%{biases.anyA{block},biases.notany{block},biases.anyB{block}};%, 
% trialblocks2use = {};
target2use = {'all'};%[1],[2],[3],[4],[5],[6],[7],[8]};%,[6]};%,[2 6],[3 7],[4 8]};%{1,2,3,4,5,6,7,8};%{1,2,3,4,5,6,7,8};
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
            biascomp = zeros(length(targnums),length(timeper));
            for tn = 1:length(targnums)
                biascomp(tn,:) = daxis(targnums(tn),timeper);
            end
            
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
