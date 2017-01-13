files2run = {'C:/Users/limblab/Desktop/datfolder/Cisek/Chewie_10262016_M.mat',...
             'C:/Users/limblab/Desktop/datfolder/Cisek/Chewie_10282016_M.mat',...
             'C:/Users/limblab/Desktop/datfolder/Cisek/Chewie_10312016_M.mat'};
         
 ABmemav = [];
 
 for fi = 1:length(files2run) 
     load(files2run{fi});
     
     mem = M.CLENS(3):M.CLENS(4);
     
     for bl = 1:length(M.pmd)
        ABmemav.all{fi,bl} = squeeze(nanmean(M.pmd{bl}.d([5 1],mem,:)./repmat(nanmean(M.pmd{bl}.d([3 7],mem,:),1),2,1,1),2))';
        correct = find(round(circ_dist(M.tt{bl}(:,15),M.tt{bl}(:,19))) == 0);
        wrong = find(round(circ_dist(M.tt{bl}(:,15),M.tt{bl}(:,19))) ~= 0);
        
        ABmemav.correct{fi,bl} = squeeze(nanmean(M.pmd{bl}.d([5 1],mem,correct)./repmat(nanmean(M.pmd{bl}.d([3 7],mem,correct),1),2,1,1),2))';
        ABmemav.wrong{fi,bl} = squeeze(nanmean(M.pmd{bl}.d([5 1],mem,wrong)./repmat(nanmean(M.pmd{bl}.d([3 7],mem,wrong),1),2,1,1),2))';
     end
     
 end
 %%
ABmem = [];
for i = 1:4
    ABmem.all{i} = cell2mat(ABmemav.all(:,i));
    ABmem.correct{i} = cell2mat(ABmemav.correct(:,i));
    ABmem.wrong{i} = cell2mat(ABmemav.wrong(:,i));
end
%%
figure; hold on; 
lims = [floor(10*min(reshape(cell2mat(ABmem.all'),[],1)))./10,ceil(10*max(reshape(cell2mat(ABmem.all'),[],1)))./10];
for i = 1:4
    subplot(1,4,i); hold on; 
    if ~isempty(ABmem.correct{i})
        plot(ABmem.correct{i}(:,1),ABmem.correct{i}(:,2),'b.');
    end
    if ~isempty(ABmem.wrong{i})
        plot(ABmem.wrong{i}(:,1),ABmem.wrong{i}(:,2),'r.');
    end
    plot(lims,[1 1],'k'); 
    plot([1 1],lims,'k');
    xlim(lims); ylim(lims);
    axis square
end