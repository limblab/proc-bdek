for j = 1:3
        
    for k = 1:length(F{1}{1});

        COMBFIRE{j}{k} = [F{1}{j}{k}; F{2}{j}{k}; F{3}{j}{k}];
    end
end

for i =1 :length(COMBFIRE)
    
    AVFIRE{i} = cell2mat(cellfun(@(x) nanmean(x,1),COMBFIRE{i},'UniformOutput',0)');
end

mincol = min(cellfun(@(x) min(x(:)),AVFIRE));
maxcol = max(cellfun(@(x) max(x(:)),AVFIRE));

figure; hold on; 
for i = 1:length(COMBFIRE)
    subplot(1,3,i);
    imagesc(AVFIRE{i},[mincol,maxcol]);
    
    for j = 1:size(AVFIRE{1},2)
       
        [~,~,~,p_bias(i,j),p_depth(i,j)]= VM_fit3(spatial_cents,AVFIRE{i}(:,j),-pi:0.01:pi);
        
    end
    
end
%%
figure; hold on; 
cs2p = {'k','b','r'};
for i = 1:length(COMBFIRE)
    plot(p_bias(i,:),cs2p{i},'LineWidth',3);
end
figure; hold on; 
cs2p = {'k','b','r'};
for i = 1:length(COMBFIRE)
    plot(p_depth(i,:),cs2p{i},'LineWidth',3);
end

%%
figure; hold on; 
for i = 1:length(COMBFIRE)
    subplot(1,3,i);
    imagesc(AVFIRE{i}-repmat(p_bias(i,:),17,1),.5*[mincol,maxcol]);
    
end
