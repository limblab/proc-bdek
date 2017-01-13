b1 = zeros(size(FA{2},2),3);
pd = zeros(size(FA{2},2),1);
BG1 = zeros(size(FA{2},2),2);
for n = 1:size(FA{2},2)
    b1(n,:) = glmfit([cos(alldays(1).tt(:,10)),sin(alldays(1).tt(:,1))],FA{1}(:,good_neurs{1}(n)),'poisson');
    
    pd(n) = atan2(b1(n,3),b1(n,2));

    BG1(n,:) = [b1(n,1) sqrt(b1(n,2)^2 + b1(n,3)^2)];
    
end
BG1g = BG1;
%%
blck = 2;
TT = alldays(blck).tt;
F = FA{blck};
conds = flipud(unique(TT(:,3)));

fitfunc = @(x) glmfit(x(:,1),x(:,2),'poisson');

[BG,BGL,BGH] = deal(cell(length(conds),1));
for i = 1:length(conds)
    
    inds = find(TT(:,3)==conds(i));
    
    for n = 1:size(F,2)     
        clc; fprintf('condition: %d/%d\nneuron %d/%d\n',i,length(conds),n,size(F,2));

        BG{i}(n,:) = glmfit(cos(TT(inds,10) - pd(n)),F(inds,n),'poisson');
        
        [BGL{i}(n,:),BGH{i}(n,:)] = boot_bounds(1000,fitfunc,[cos(TT(inds,10)-pd(n)) F(inds,n)],2.5,97.5);
    end

end
%%
[badinds,dBase,dGain] = deal(cell(length(conds),1));
BGnew = BG;
BGLnew = BGL;
BGHnew = BGH;
for i = 1:length(conds)
    %badinds{i} = find(sum(abs(BG{i}),2) > 10e4);
    badinds{i} = find(sum(F)<150);
    BGnew{i}(badinds{i},:) = nan;
    BGLnew{i}(badinds{i},:) = nan;
    BGHnew{i}(badinds{i},:) = nan;
    
    MINS_U{i} = exp(BGnew{i}(:,1) - BGnew{i}(:,2));
    MAXS_U{i} = exp(BGnew{i}(:,1) + BGnew{i}(:,2));
    
end
MINS_C = exp(BG1g(:,1) - BG1g(:,2));
MAXS_C = exp(BG1g(:,1) + BG1g(:,2));

[dB,dG] = deal(zeros(length(conds),2));
figure; hold on; cols = {'b','r'};
%plot(MINS_C,MAXS_C,'k.');
for i = 1:length(conds)
    
%     dBase{i} = (exp(BGnew{i}(:,1)-BG1g(:,1)))./exp(BG1g(:,1));
%     dGain{i} = (exp(BGnew{i}(:,2)-BG1g(:,2)))./exp(BG1g(:,2));
%     
%     dB(i,2) = nanmean(dBase{i});
%     dG(i,2) = nanmean(dGain{i});
%     [dB(i,1),dB(i,3)] = boot_bounds(1000,@nanmean,dBase{i},2.5,97.5);
%     [dG(i,1),dG(i,3)] = boot_bounds(1000,@nanmean,dGain{i},2.5,97.5);
    
    plot((MINS_U{i}-MINS_C)./MINS_C,(MAXS_U{i}-MAXS_C)./MAXS_C,'.','Color',cols{i});
    plot(nanmean((MINS_U{i}-MINS_C)./MINS_C),nanmean((MAXS_U{i}-MAXS_C)./MAXS_C),'o','Color',cols{i},'MarkerSize',14);
    
    %plot((BGnew{i}(:,1)-BG1g(:,1))./BG1g(:,1),(BGnew{i}(:,2)-BG1g(:,2))./BG1g(:,2),'.','Color',cols{i});
%     for n = 1:size(BG1g,1)
%         %plot(BGnew{i}(:,1)-BG1(:,1),BGnew{i}(:,2)-BG1(:,2),'.','Color',cols{i})
%         basediff_1 = [BGLnew{i}(n,1)-BG1g(n,1) BGHnew{i}(n,1)-BG1g(n,1)];
%         gaindiff_1 = [BGnew{i}(n,2)-BG1g(n,2) BGnew{i}(n,2)-BG1g(n,2)];
%         
%         basediff_2 = [BGnew{i}(n,1)-BG1g(n,1) BGnew{i}(n,1)-BG1g(n,1)];
%         gaindiff_2 = [BGLnew{i}(n,2)-BG1g(n,2) BGHnew{i}(n,2)-BG1g(n,2)];
%         
%         
%         plot(basediff_1,gaindiff_1,'-','Color',cols{i});
%         plot(basediff_2,gaindiff_2,'-','Color',cols{i});
% %     
% 
%     end
    
%     plot([dB(i,1) dB(i,3)],[dG(i,2),dG(i,2)],'-','Color',cols{i});
%     plot([dB(i,2) dB(i,2)],[dG(i,1),dG(i,3)],'-','Color',cols{i});
end

% for n = 1:size(BG1g,1)
%     
%     plot([MINS_C(n) cellfun(@(x) x(n), MINS_U)],[MAXS_C(n) cellfun(@(x) x(n), MAXS_U)],'k');
% end

