function [pval,pval_o,okinds,percstd] = match_dists(RT,kstd,varargin)

RT(RT < 0) = nan;

Qmu = nanmean(RT(:,1));
if nargin>1
    perclist = 1;
else
    kstd = std(RT(:,1));
    perclist = 1:-0.05:0.2;
end

liklist = unique(RT(:,2));
[ds,newds,keepers,Globinds,keepinds] = deal(cell(length(liklist),1));

stoploop = 0;
for pl = 1:length(perclist)
    
    if stoploop == 0
        Qstd = perclist(pl).*kstd;

        for i = 1:length(liklist)
            Globinds{i} = find(RT(:,2)==liklist(i));
            ds{i} = RT(Globinds{i},1);

            P = pdf('normal',ds{i},nanmean(ds{i}),nanstd(ds{i}));
            Q = pdf('normal',ds{i},Qmu,Qstd);

            keepers{i} = (Q./P > 1) | (Q./P > rand(length(P),1));

            newds{i} = ds{i}(keepers{i});
            keepinds{i} = Globinds{i}(keepers{i});
        end
        [diffmeans,pval] = ttest2(newds{1},newds{2});
        [~,pval_o] = ttest2(ds{1},ds{2});

        if diffmeans == 0
            stoploop = 1;
            percstd = perclist(pl);
        end
    end
end
    
    
okinds = sortrows([keepinds{1}; keepinds{2}]);
