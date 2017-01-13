ABt = cell(4,1);

tis = cell(4,8);
for i = 1:4
    for j = 1:8 
        tis{i,j} = find(reachids{i}==j);
    end
    
    for j = 1:4
        ABt{j,i} = [reshape(M.pmd{i}.d([5 1],mem,tis{i,j})./repmat(nanmean(M.pmd{i}.d([3 7],mem,tis{i,j}),1),2,1,1),2,[])';...
                    reshape(M.pmd{i}.d([1 5],mem,tis{i,j+4})./repmat(nanmean(M.pmd{i}.d([3 7],mem,tis{i,j+4}),1),2,1,1),2,[])'];
    end
%     ABt{:,i} = cell2mat(ABt());
end
%%
all = reshape(cell2mat(reshape(ABt,[],1)),[],1);
lims = [floor(10*min(all))./10,ceil(10*max(all))./10];
colrs = {'b','r','g','m'};
figure; hold on; 
for i = 1:4; 
    subplot(1,4,i); hold on; 
    for j = 1:4
        plot(ABt{j,i}(:,1),ABt{j,i}(:,2),'.','Color',colrs{j}); 
    end
    plot(lims,[1 1],'k'); 
    plot([1 1],lims,'k');
    xlim(lims); ylim(lims);
    axis square
end