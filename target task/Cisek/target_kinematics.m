POS = bdfM.pos; POS(:,2) = POS(:,2)-3; POS(:,3) = POS(:,3)+33;
conds = unique(tt_pred(:,20));
projec = cell(length(conds),1);
for lik = 1:length(conds);
    
    trials = find(tt_pred(:,20)==conds(lik));
    
    for i = 1:length(trials)
        
        t1 = find(POS(:,1) > tt_pred(trials(i),9),1,'first');
        t2 = find(POS(:,1) > tt_pred(trials(i),9)+2,1,'first');
        
        ps = POS(t1:t2,2:3);
        
        movevec = [cos(tt_pred(trials(i),17)) sin(tt_pred(trials(i),17))];    
        moverep = repmat(movevec,size(ps,1),1);
 
        dot1 = dot(ps,moverep,2); 
        
        projec{lik}(i,:) = dot1';
        clc; fprintf('Condition: %d/%d\nTrial: %d/%d\n',lik,length(conds),i,length(trials));
       
    end
end
%%
colors2plot = distinguishable_colors(length(projec));

Mprojs = cellfun(@mean,projec,'UniformOutput',0);
PROJ = (vertcat(Mprojs{:}))';

figure; hold on;
for i= 1:length(projec)
    plot(PROJ(:,i),'Color',colors2plot(i,:));
end

%%
[projecALL,deltaproj] = deal(zeros(size(tt_pred,1),2001));
RT = zeros(size(tt_pred,1),1);
for i = 1:size(tt_pred,1)
        
        t1 = find(POS(:,1) > tt_pred(i,9),1,'first');
        t2 = find(POS(:,1) > tt_pred(i,9)+2,1,'first');
        
        ps = POS(t1:t2,2:3);
        
        movevec = [cos(tt_pred(i,17)) sin(tt_pred(i,17))];    
        moverep = repmat(movevec,size(ps,1),1);
 
        dot1 = dot(ps,moverep,2); 
        
        projecALL(i,:) = dot1';
        deltaproj(i,:) = [0 diff(dot1')];
        
        takeoff = find(abs(diff(dot1)) > 0.001,1,'first');
        if isempty(takeoff); takeoff = 1; end
        RT(i) = takeoff;
        
        clc; fprintf('Trial: %d/%d\n',i,size(tt_pred,1));
end