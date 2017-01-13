
prob_at_end = chi2cdf(4*n*(tt1^2 - n^2*r^2)/(tt1^2 - 2*n^2),1);
prob_at_end = chi2cdf(5*log((1-(r*cos(tt))^2)/(1-r^2)),1);
%%
[mid_block,low_block,high_block] = deal(zeros(size(midbound{1},1),length(midbound)));
for i = 1:length(midbound)
    
    mid_block(:,i) = midbound{i}(:,end);
    low_block(:,i) = lowerbound{i}(:,end);
    high_block(:,i) = upperbound{i}(:,end);
end

blockx = 1:size(mid_block,2);
figure; hold on;
for i = 1:size(mid_block,1)
    
    plot(blockx,mid_block(i,:),colors_of_plot{i});
    patch([blockx fliplr(blockx)],[low_block(i,:) fliplr(high_block(i,:))],colors_of_plot{i},'FaceAlpha',.5,'EdgeAlpha',.5);
end