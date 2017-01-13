trialkin = cell(length(alldays),1);
tseq = alldays(1).kin.pos(:,1);
rdir = cell(length(alldays),1);
todir = cell(length(alldays),1);
toang = cell(length(alldays),1);

for i = 1:length(alldays)
    
    for j = 1:size(TT{i},1)
        
        st = TT{i}(j,6);
        et = TT{i}(j,3);
        
        trialkin{i}{j,:} = alldays(1).kin.pos(tseq>st & tseq<et,2:3);
        if et>st
            rdir{i}(j,:) = atan2(trialkin{i}{j}(end,2),trialkin{i}{j}(end,1));
        else
            rdir{i}(j,:) = NaN;
        end
        
        tot = TT{i}(j,12);
        totrace = alldays(1).kin.pos(tseq>tot & tseq<(tot+.1),2:3);
        if ~isnan(tot)
            to = atan2(totrace(end,2)-totrace(1,2),totrace(end,1)-totrace(1,1));
            cor = atan2(trialkin{i}{j}(end,2)-totrace(1,2),trialkin{i}{j}(end,1)-totrace(1,1));
            todir{i}(j,:) = circ_dist(to,cor); 
            toang{i}(j,:) = to;
        else
            todir{i}(j,:) = NaN;
            toang{i}(j,:) = NaN;
        end
        
    end
end
todirA = cell(length(rdir),1);
toangA = cell(length(rdir),1);
p = polyfit(rdir{1},todir{1},8);
for i = 1:length(todir)
    for j = 1:length(todir{i})
        todirA{i}(j,:) = circ_dist(todir{i}(j),polyval(p,rdir{i}(j))); 
        toangA{i}(j,:) = circ_dist(toang{i}(j),polyval(p,rdir{i}(j)));
    end
end

%%
figure; hold on; 

for i = 1:length(trialkin{1})
    plot(trialkin{1}{i}(:,1),trialkin{1}{i}(:,2),'b');
end

for i = 1:length(trialkin{2})
    plot(trialkin{2}{i}(:,1),trialkin{2}{i}(:,2),'r');
end
