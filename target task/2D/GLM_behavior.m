endpoints = alldays(2).tt(:,10);

% Likelihood basis functions 
thets = -pi:0.001:pi;
basisK = 2.^(0:1:9);

b = zeros(length(basisK),length(thets));
for nbasis = 1:length(basisK)
    bK = basisK(nbasis);
    b(nbasis,:) = circ_vmpdf(thets,0,bK);
end
b = (b - repmat(min(b,[],2),1,size(b,2)))./repmat(max(b,[],2) - min(b,[],2),1,size(b,2));

% Prior basis functions
tfilt = 1./[10000 1000 100 10 1];

d = zeros(length(tfilt),10000);
for nbasis = 1:length(tfilt)
    d(nbasis,:) = exp(-tfilt(nbasis).*(1:1:10000));
end

thetas = 0:0.001:2*pi;
MLE = zeros(length(alldays(2).tt),length(basisK)*length(tfilt));
MODEL = zeros(length(basisK)*length(tfilt),2);
MLE(i,:) = NaN;
for i = 2:length(alldays(2).tt)
    
    slicelist = mod(alldays(2).slices(i,:),2*pi);
    
    slicetrain = zeros(1,length(thets));
    slicetrain(round(1000*slicelist)) = 1;
    repslice = repmat(slicetrain,1,2);
    %repslice = [zeros(1,2*length(thets)) slicetrain zeros(1,2*length(thets))];
    
    likeli_est = zeros(length(basisK),length(repslice)./2);
    for j = 1:length(basisK)
        likeli_long = conv(repslice,b(j,:),'same'); 
        
        likeli_est = likeli_long(1:end/2);
        likeli_est = (likeli_est - min(likeli_est))/...
            (max(likeli_est) - min(likeli_est));
        
        prev_trials = alldays(2).tt(1:i-1,10)';   

        prior_est = zeros(length(tfilt),length(thets));
        for k = 1:length(tfilt)
            scaled_prev = repmat(d(k,1:length(prev_trials)),2,1).*...
                [cos(fliplr(prev_trials)); sin(fliplr(prev_trials))];
            prior_sum = sum(scaled_prev,2);
            prior_est = .5.*length(prev_trials)*(1 + ...
                norm(prior_sum).*cos(thetas - atan2(prior_sum(2),prior_sum(1)))); 

            post_est = likeli_est.*prior_est;
            
            MLE(i,(j-1)*length(tfilt) + k) = thetas(post_est == max(post_est));
            if i == 2; MODEL((j-1)*length(tfilt) + k, :) = [basisK(j) tfilt(k)]; end 
            
        end
        
    end

end
% % %
% [betas, dev, stats] = glmfit(MLE(2:end,:),alldays(2).tt(2:end,10));
% yhat = glmval(betas,MLE,'identity'); 

%%
[P,S,yhat] = deal(cell(size(MLE,2),1));
[SSR,SST,R2] = deal(zeros(size(MLE,2),1));
for m = 1:size(MLE,2)
    [P{m},S{m}] = polyfit(MLE(2:end,m),endpoints(2:end),1);
    
    yhat{m} = P{m}(2) + P{m}(1)*MLE(2:end,m);
   
    SSR(m) = sum((yhat{m} - endpoints(2:end)).^2);
    SST(m) = sum((mean(endpoints(2:end)) - endpoints(2:end)).^2);
    R2(m) = 1 - SSR(m)/SST(m);
end
