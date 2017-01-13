figure; hold on; 
cp = {'k','b','r'};
[kapps,confid,cmeans,cstds] = deal(cell(length(BOOTDIR),1));
[total_kapps,total_means,total_conf] = deal(zeros(length(alldays(2).tt),size(BOOTDIR{1}{1},2)));
for i = 1:length(BOOTDIR)
    
    for j = 1:length(BOOTDIR{i})
%         
%         patch([1:size(BOOTDIR{i}{j},2) size(BOOTDIR{i}{j},2):-1:1],...
%               [BOOTDIR{i}{j}(1,:) fliplr(BOOTDIR{i}{j}(2,:))],cp{i},'FaceAlpha',0.5,'EdgeColor','none');
%   
        for tb = 1:size(BOOTALL{i}{j},2)
            kapps{i}(j,tb) = circ_kappa(BOOTALL{i}{j}(:,tb));
            cmeans{i}(j,tb) = circ_mean(BOOTALL{i}{j}(:,tb));
            cstds{i}(j,tb) = circ_std(BOOTALL{i}{j}(:,tb));
        end
        confid{i}(j,:) = diff(BOOTDIR{i}{j});
          
    end
    
%     plot(nanmean(confid{i}),cp{i}); 
%     [bl,bh] = boot_bounds(1000,@nanmean,confid{i},2.5,97.5);
    
%     plot(nanmean(log(kapps{i})),cp{i}); 
%     [bl,bh] = boot_bounds(1000,@nanmean,log(kapps{i}),2.5,97.5);

    plot(circ_mean(abs(cmeans{i})),cp{i}); 
    [bl,bh] = boot_bounds(1000,@circ_mean,abs(cmeans{i}),2.5,97.5);

%       plot(nanmean(cstds{i}),cp{i}); 
%       [bl,bh] = boot_bounds(1000,@nanmean,cstds{i},2.5,97.5);

    
    patch([1:size(BOOTDIR{i}{j},2) size(BOOTDIR{i}{j},2):-1:1],...
              [bl' fliplr(bh')],cp{i},'FaceAlpha',0.5,'EdgeColor','none');

    total_kapps(like_ind{1}{i},:) = kapps{i};
    total_means(like_ind{1}{i},:) = cmeans{i};
%     total_conf(like_ind{1}{i},:) = nanmean(confid{i},2);

end
%%
figure; hold on;
c = colormap; 
for i = 1:length(alldays(2).tt)
    
    plot(alldays(2).tt(i,9),alldays(2).tt(i,10),'.','Color',c(width2col(total_conf(i)),:),'MarkerSize',5);
end


%%
[bootkap,bootvar] = deal(zeros(size(alldays(2).slices,1),1));
for i = 1:size(alldays(2).slices,1)
    
    [~,~,~,bdist] = boot_bounds(1000,@circ_mean,alldays(2).slices(i,:),2.5,97.5);
    bootkap(i,:) = circ_kappa(bdist);
end

%%
thets = -pi:0.01:pi;
for i = 1:size(alldays(2).tt,1)
    
    pridist = circ_vmpdf(thets,alldays(2).tt(i,9),total_conf(i));
end
