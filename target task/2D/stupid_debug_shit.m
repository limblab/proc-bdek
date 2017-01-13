[PON,DS,H,CI] = deal(cell(3,1));
for i = 1:length(VA)
    for j = 1:length(VA{i}) 
        for j2 = 1:(length(VA{i}{j})-1)
            k = j2 +1;
            
            d1 = VA{i}{j}{k}(LI{i}{1},:);
            d2 = VA{i}{j}{k}(LI{i}{2},:);

            nd1 = sum(~isnan(d1));
            nd2 = sum(~isnan(d2));

            for u = 1:size(d1,2);
%                 [H{j2}{j}{i}(u,:),~,CI{j2}{j}{i}(u,:)] = ttest2(d2(:,u),d1(:,u),0.05,'both','unequal');

                if nd1(u)>0 && nd2(u)>0
                    [~,hreturn] = ranksum(d2(:,u),d1(:,u));
                    H{j2}{j}{i}(u,:) = double(hreturn);
                else
                    H{j2}{j}{i}(u,:) = NaN;
                end
            end
            validunits = ~isnan(H{j2}{j}{i});
            nvalidunits = sum(validunits);
            
            incs = nanmean(d2)>nanmean(d1);
            decs = nanmean(d2)<nanmean(d1);
                       
            pover = nansum(H{j2}{j}{i}==1 & incs'==1);%./nvalidunits;
            punder = nansum(H{j2}{j}{i}==1 & decs'==1);%./nvalidunits;
            pnon = sum(nd1>5 & nd2>5);%nansum(H{j2}{j}{i}==0)./nvalidunits;
 
            PON{j2}{j}(i,:) = [pover punder pnon];
                        
        end
    end
end