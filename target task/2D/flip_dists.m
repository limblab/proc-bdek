function [okinds,pval,pval_o,percstd] = flip_dists(RT)

RT(RT < 0) = nan;

% if nargin>1
%     perclist = 1;
% else
%     kstd = std(RT(:,1));
%     perclist = 1:-0.05:0.2;
% end
perclist = 1:-0.05:0.2;

liklist = unique(RT(:,2));
[Qmu,Qstd] = deal(zeros(length(liklist),1));
for i = 1:length(liklist)
    Qmu(i) = nanmean(RT(RT(:,2)==liklist(i),1));
    Qstd(i) = nanstd(RT(RT(:,2)==liklist(i),1));
end
Qmu = flipud(Qmu);
Qstd = flipud(Qstd);

[ds,newds,keepers,Globinds,keepinds] = deal(cell(length(liklist),1));
stoploop = 0;
for pl = 1:length(perclist)
    
    if stoploop == 0

        for i = 1:length(liklist)
            Globinds{i} = find(RT(:,2)==liklist(i));
            ds{i} = RT(Globinds{i},1);

            P = pdf('normal',ds{i},nanmean(ds{i}),nanstd(ds{i}));
            Q = pdf('normal',ds{i},Qmu(i),Qstd(i)*perclist(pl));

            keepers{i} = (Q./P > 1) | (Q./P > rand(length(P),1));

            newds{i} = ds{i}(keepers{i});
            keepinds{i} = Globinds{i}(keepers{i});
        end
        [diffmeans,pval] = ttest2(newds{1},newds{2});
        [difforiginal,pval_o] = ttest2(ds{1},ds{2});
        
        signdif = sign(mean(newds{1})-mean(newds{2}));
        signdif_O = sign(mean(ds{1})-mean(ds{2}));

        if (diffmeans == difforiginal) && (signdif~=signdif_O);
            stoploop = 1;
            percstd = perclist(pl);
        end
    end
end
    
    
okinds = sortrows([keepinds{1}; keepinds{2}]);
