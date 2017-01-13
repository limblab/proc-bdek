b_units = alldays(1).PMd_units;

tend = max(alldays(2).tt(:,10));
contFRbin = zeros(length(b_units),ceil(tend*1000/50)-1);
for i = 1:length(b_units)
    
    clc; fprintf('%d/%d\n',i,length(b_units));
    [rast,t,idx] = histcounts(b_units{i}(2:end),0:0.05:tend);
    t2 = mean([t(1:end-1); t(2:end)]);
    contFR = train2cont(rast,2)/50;
    
    contFRbin(i,:) = contFR;%bin_array(contFR,1,round((length(contFR)+1)/50),'mean');

end
%%
binTT = alldays(2).tt;
binTT(:,[4:10 20]) = round(binTT(:,[4:10 20])/0.05);

%%
[pcs,wghts] = pca(contFRbin');

%%
trl = cell(size(binTT,1),1);

for i = 1:size(binTT,1)
    
    trl{i} = wghts(binTT(i,6):binTT(i,7),:)';
end

%%
figure; hold on;
clrs = distinguishable_colors(8);
for i = 1:length(t1)
    plot(trl{i}(1,:),'Color',clrs(test_targs(t1(i)),:));
end