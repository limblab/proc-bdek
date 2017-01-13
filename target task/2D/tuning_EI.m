%% define as excitatory or inhibitory
PrefD = cell(length(neurons{1}),1);
inhibit = []; excite = [];
ex_in = zeros(length(neurons{1}),1);
for i = 1:length(neurons{1})
    
    AV1 = mean(neurons{1}{i}.tuning);
    AV2 = mean(neurons{2}{i}.tuning);
    
    MIN2 = min(neurons{2}{i}.tuning);
    MAX2 = max(neurons{2}{i}.tuning);
    
%     if (AV2-AV1)>0
%         PrefD{i} = wrapped_cents(neurons{1}{i}.tuning==max(neurons{2}{i}.tuning));
%         %excite = [excite;i];
%     elseif (AV2-AV1)<0 
%         PrefD{i} = wrapped_cents(neurons{1}{3}{i}.tuning==min(neurons{2}{i}.tuning));
%         %inhibit = [inhibit;i];
%     else 
%         PrefD{i} = NaN;
%         
%     end
%     
%     if length(PrefD{i}) > 1
%         PrefD{i} = NaN;
%     end
    
    PDind = find(abs(neurons{2}{i}.tuning-AV1)==max(abs(neurons{2}{i}.tuning-AV1)));
    PrefD{i} = wrapped_cents(PDind);
    
    if length(PrefD{i}) > 1
        PrefD{i} = NaN;
    end
    
    if neurons{2}{i}.tuning(PDind(1))-AV1 > 0 && ~isnan(PrefD{i})
        ex_in(i) = 1;
    elseif neurons{2}{i}.tuning(PDind(1))-AV1 < 0 && ~isnan(PrefD{i})
        ex_in(i) = -1;
    else
        ex_in(i) = NaN;
    end
    
    if MAX2<AV1
        inhibit = [inhibit;i];
    elseif MIN2>AV1
        excite = [excite;i];
    end
    
end
PrefDs = vertcat(PrefD{:});

nears = find(abs(circ_dist(PrefDs,circ_mean(alldays(2).tt(:,10))))<pi/2);
fars  = find(abs(circ_dist(PrefDs,circ_mean(alldays(2).tt(:,10))))>pi/2);

far_in = fars(ismember(fars,find(ex_in==-1)));
far_ex = fars(ismember(fars,find(ex_in==1)));
near_in = nears(ismember(nears,find(ex_in==-1)));
near_ex = nears(ismember(nears,find(ex_in==1)));