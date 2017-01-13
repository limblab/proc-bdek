smrast1 = smrast(1:length(alldays(1).tt));
smrast2 = smrast((length(alldays(1).tt)+1):end);

%% separate by neuron
neurarray = cell(length(array),1);
[loo_tbt,loo_ax,loo_scores] = deal(cell(length(smrast1),1));
for i = 1:length(smrast1)
    clc; fprintf('%d/%d\n',i,length(smrast1));
    ileftout = smrast1; ileftout(i) = [];
    for j = 1:length(array)
         neurarray{j} = cell2mat(transpose(cellfun(@(x) x(j,:),ileftout,'UniformOutput',0)));
    end
    
    totarray = vertcat(neurarray{:})';
    [loo_ax{i}, loo_scores{i}] = pca(totarray);
    
    loo_proj = smrast1{i}'*loo_ax{i};
    loo_proj(:,11:end) = [];
    
    loo_tbt{i} = loo_proj;
end

%% Trim traces to relevant (and equal length) region
kpinds = 1:(750/dsamplesize+1);
LOOtrimmed = cell(10,1);
for i = 1:length(loo_tbt) % trial
    for j = 1:10 % PC
        LOOtrimmed{j}(i,:) = loo_tbt{i}(kpinds,j)';
    end
end
%%

[loo_b] = deal(cell(10,1));
for j = 1:10
    for i = 1:size(LOOtrimmed{1},1)
        Xctbt = [BLCK1{j}(targid1(i):end,:) ; BLCK1{j}(1:(targid1(i)-1),:)]';
        loo_b{j}(i,:) = Xctbt\(LOOtrimmed{j}(i,:)');
    end
end