function [preds,posts,truelabs] = naive_Bayes_cv(data,labels,numcv)

if isempty(numcv) % Do leave-one-out
    numcv = length(labels);
end
%% Shuffle indices for data and labels
N = length(labels); % length of data
shuff_inds = randperm(N)'; % shuffle trials

%% Partition for cross validation
cv_inds = round(linspace(1,(N+1),numcv+1)); 

%% Loop through each fold of cross-validation
[testi,traini,posterior,predict,truelabels] = deal(cell(numcv,1));
for i = 1:numcv
    clc; fprintf('Fold %d/%d\n',i,numcv);
    
    testi{i} = shuff_inds(cv_inds(i):(cv_inds(i+1)-1));
    traini{i} = shuff_inds(~ismember(shuff_inds,testi{i}));
    
    % Fit classifier
    nb = NaiveBayes.fit(data(traini{i},:),labels(traini{i}));
    
    % Test on test data
    posterior{i} = nb.posterior(data(testi{i},:));
    predict{i} = nb.predict(data(testi{i},:));
    truelabels{i} = labels(testi{i});
   
end
    
posts_shuff = vertcat(posterior{:});
preds_shuff = vertcat(predict{:});
truelabs_shuff = vertcat(truelabels{:});

posts(shuff_inds,:) = posts_shuff;  
preds(shuff_inds,:) = preds_shuff;
truelabs(shuff_inds,:) = truelabs_shuff; 
    
    
