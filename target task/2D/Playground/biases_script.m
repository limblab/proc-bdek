region = pretarg;
[biases] = deal(cell(2,1));
biases_full = cell(2,8);
for i = 1:8
    biases{1}(i,:) = nanmean(diff(M.pmd{2}.d([1 5],region,trials.corrects{2}.targs{i}),[],1),3); 
    biases{2}(i,:) = nanmean(diff(M.m1{2}.d([1 5],region,trials.corrects{2}.targs{i}),[],1),3);
    
    biases_full{1,i} = squeeze(diff(M.pmd{2}.d([1 5],region,trials.corrects{2}.targs{i}),[],1)); 
    biases_full{2,i} = squeeze(diff(M.m1{2}.d([1 5],region,trials.corrects{2}.targs{i}),[],1));
end

figure; hold on; 
symbl = {'.','o','^','+'};

[cor,r2,av,avf] = deal(zeros(4,2));
for i = 1:4
    plot(biases{1}(i,:),biases{2}(i,:),symbl{i},'Color','b','MarkerSize',10);
    plot(biases{1}(i+4,:),biases{2}(i+4,:),symbl{i},'Color','r','MarkerSize',10);

    b1 = reshape(biases{1}([i i+4],:),[],1);
    b2 = reshape(biases{2}([i i+4],:),[],1);
    
    b1f = reshape(biases_full{1,i},[],1);
    b2f = reshape(biases_full{2,i},[],1);
    nanlocs = isnan(b1f+b2f);
    b1f(nanlocs) = [];
    b2f(nanlocs) = [];
    
    cc = corrcoef(b1,b2);
    ccf = corrcoef(b1f,b2f);
    
    cor(i,1) = cc(1,2);
    cor(i,2) = ccf(1,2);
    
    r2(i,1) = cor(i,1).^2;
    r2(i,2) = cor(i,2).^2;
    
    av(i,:) = [nanmean(b1),nanmean(b2)];
    avf(i,:) = [nanmean(b1f),nanmean(b2f)];


end


