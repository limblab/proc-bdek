ABmemtime = cell(1,4);
for bli = 1:4
    for i = 1:size(M.pmd{2}.d,2)
        ABmemtime{bli}(:,:,i) = reshape(M.pmd{bli}.d([5 1],i,:)./repmat(nanmean(M.pmd{bli}.d([3 7],i,:),1),2,1,1),2,[])';
    end
end
%%
% Remove bad times
removetimes = 67:79;
xcourse = 1:size(ABmemtime{1},3);
for i = 1:length(ABmemtime)
    ABmemtime{i}(:,:,removetimes) = NaN;
end
xcourse(removetimes) = NaN;

%% Add sister plot
Y(removetimes,:) = NaN;


%%
block = 4;


check_epoch = @(x) find(M.CLENS<=x,1,'last');
figure; hold on;
plot([0.5 3],[1 1],'k'); 
plot([1 1],[0.5 3],'k'); 
xlim([.5 3]); ylim([.5 3]); axis square;

colrs = {'k','b','b','r','r','g','g','g'};
correct = find(round(circ_dist(M.tt{block}(:,15),M.tt{block}(:,19))) == 0);
wrong = find(round(circ_dist(M.tt{block}(:,15),M.tt{block}(:,19))) ~= 0);
all = 1:size(M.tt{block},1);

trialset = all; 

figure; hold on; pAB = subplot(1,2,2); hold on; 
plot(pAB,[0.5 3],[1 1],'k'); 
plot(pAB,[1 1],[0.5 3],'k'); 
xlim([.5 3]); ylim([.5 3]); axis square;

pT = subplot(1,2,1); hold on;
plot(pT,[0 M.CLENS(end)],[1 1],'k');
xlim([0 M.CLENS(end)]); 
ylim([0.7 1.7]); 

for i = 1:size(ABmemtime{block},3)
    
    if ~isnan(xcourse(i))
        cur_ep = check_epoch(i);
    %     a = plot(squeeze(ABmemtime{block}(:,1,i)),squeeze(ABmemtime{block}(:,2,i)),'.','Color',colrs{cur_ep}); 
        xs = reshape(ABmemtime{block}(trialset,1,M.CLENS(cur_ep):i),1,[]);
        ys = reshape(ABmemtime{block}(trialset,2,M.CLENS(cur_ep):i),1,[]);
        a = plot(pAB,xs,ys,'.','Color',colrs{cur_ep}); 
        
        
        b = plot(pT,Y(1:i,1),'b-','LineWidth',5); 
        c = plot(pT,Y(1:i,2),'r-','LineWidth',5);
        
        
        if ismember(i,M.CLENS(2:(end-2))-1)
            plot(pT,(i+1)*[1 1],[0.7 1.7],'k'); 
            pause;
        else
            pause(0.1);
        end
        if i~= size(ABmemtime{block},3)
            delete(a);
        end
%         delete(b); delete(c);
    end 
end