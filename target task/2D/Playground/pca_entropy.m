bszs = [0.1:0.1:20];

for q = 1:length(bszs)
ensemble = alldays(1).M1_units;
binsize = bszs(q);
T = alldays(1).tt(end,7);
N = length(ensemble);

timebins = 1:1/binsize:T;
rates = zeros(length(ensemble),length(timebins));
for n = 1:N
    
    rates(n,:) = histc(ensemble{n}',timebins)./binsize;
    
    
end
%%
dfunc = @(x) (x(3:end)-x(1:end-2))./2;

[~,weights,score] = pca(rates');

%%
vaf = sum(repmat(score,1,length(score)).*triu(ones(length(score))),1)./sum(score);

numpcs = 10;

threshold = mean(score(numpcs:end));

truncated = weights(:,1:numpcs)';

discretized = abs(truncated) > threshold; 

state = zeros(1,size(discretized,2));
for i = 1:size(discretized,2) 
    state(i) = bin2dec(num2str(discretized(:,i)'));
end
%%
pstate = zeros(1,2^numpcs);
for i = 1:2^numpcs;
    pstate(i) = sum(state==i)./length(state);
end
S_entropy(q) = -nansum(pstate.*log2(pstate));
end
