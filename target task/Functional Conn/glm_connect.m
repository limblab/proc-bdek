BINSIZE = 10;

[iht, ihbas, ihbasis] = makeRaisedCosBasis(5, BINSIZE/1000, [2/1000 25/1000]*7,1*1.5*1e-2,1);
%%
trange = alldays(1).tt(end,1);

ts = 0:0.01:trange;

s_PMd = cell(length(alldays(1).PMd_units),1);
s_M1 = cell(length(alldays(1).M1_units),1);
for i =1 :length(alldays(1).PMd_units)
    s_PMd{i} = train2bins(alldays(1).PMd_units{i},ts);
    clc; fprintf('Area: 1/2\nUnit: %d/%d\n',i,length(alldays(1).PMd_units));
end
for i = 1:length(alldays(1).M1_units)
    s_M1{i} = train2bins(alldays(1).M1_units{i},ts);
    clc; fprintf('Area: 2/2\nUnit: %d/%d\n',i,length(alldays(1).M1_units));
end

%% Compile cell of spike trains convolved with basis functions
bPMd = cell(length(s_PMd),1);
bM1 = cell(length(s_M1),1);
for j = 1:size(ihbasis,2)
    for i = 1:length(s_PMd)
        c = conv(ihbasis(:,j)',s_PMd{i});
        bPMd{i}(j,:) = c(1:length(s_PMd{i}));
        clc; fprintf('basis: %d/%d\nArea 1/2\nUnit %d/%d\n',j,size(ihbasis,2),i,length(s_PMd));
    end
    for i = 1:length(s_M1)
        c = conv(ihbasis(:,j)',s_M1{i});
        bM1{i}(j,:) = c(1:length(s_M1{i}));
        clc; fprintf('basis: %d/%d\nArea 2/2\nUnit %d/%d\n',j,size(ihbasis,2),i,length(s_M1));
    end
end
%%
allunits = [s_PMd(keepPMd); s_M1(keepM1)];
allunits_lagged = [bPMd(keepPMd); bM1(keepM1)];

% allunits = [s_PMd(pmdnarrow); s_PMd(pmdbroad)];
% allunits_lagged = [bPMd(pmdnarrow); bPMd(pmdbroad)];

predictmat = vertcat(allunits_lagged{:})';
kernels = cell(length(allunits),1);
for m = 1:length(allunits)
    [bs] = glmfit(predictmat,allunits{m}','poisson');
    %
    kernel_ws = reshape(bs(2:end),size(ihbasis,2),[]);
    for i = 1:size(kernel_ws,2)
        kernels{m}(i,:) = sum(repmat(kernel_ws(:,i),1,size(ihbasis,1)).*ihbasis');
    end
    clc; fprintf('Predicting unit %d/%d...\n',m,length(allunits));
end
%%
kernel_sum = cell2mat(cellfun(@(x) mean(x,2)',kernels,'UniformOutput',0));
optimal_lag = cell2mat(cellfun(@(x) cell2mat(cellfun(@(y) find(abs(y)==max(abs(y)),1,'first'),...
    mat2cell(x,ones(length(allunits),1),size(ihbasis,1)),'UniformOutput',0))',kernels,'UniformOutput',0));
offdiags = kernel_sum.*(1+diag(-1*ones(size(kernel_sum,1),1)));
figure;
imagesc(offdiags); hold on;
plot([0 length(allunits)],[1 1]*length(allunits)/2,'w','LineWidth',5);
plot([1 1]*length(allunits)/2,[0 length(allunits)],'w','LineWidth',5);



