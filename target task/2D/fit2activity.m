SC = trial_counts{1};
SCU = zeros(size(SC{1},1),size(SC{1},2),length(SC));
for i =1:length(SC)
    SCU(:,:,i) = SC{i}/10;
end

%%
SCC = zeros(size(TA{1},2),size(TA{1},1),length(TA));
for i = 1:length(TA)
    SCC(:,:,i) = TA{i}';
end

mint = squeeze(min(min(SCC,[],1),[],3));
maxt = squeeze(max(max(SCC,[],1),[],3));


SCU = (SCU - repmat(mint,size(SCU).^[1 0 1]))./repmat(maxt-mint,size(SCU).^[1 0 1]);
SCC = (SCC - repmat(mint,size(SCC).^[1 0 1]))./repmat(maxt-mint,size(SCC).^[1 0 1]);
%%
for i = 1:size(SCU,1)
    
    SCUT = squeeze(SCU(i,:,:));
    
    for j = 1:size(SCUT,2)
        
        SCUTt = SCUT(:,j);
        
        W * SCC