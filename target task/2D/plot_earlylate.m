hinds = find(alldays(2).tt(:,3)==min(alldays(2).tt(:,3)));

[FDE,FDL] = deal(zeros(size(firing_diffs{1}{1},2),2));

[FDE(:,1), FDE(:,2)] = boot_bounds(1000,@nanmean,firing_diffs{1}{1}(hinds,:),2.5,97.5);
[FDL(:,1), FDL(:,2)] = boot_bounds(1000,@nanmean,firing_diffs{1}{2}(hinds,:),2.5,97.5);

FDEA = mean(firing_diffs{1}{1}(hinds,:));
FDLA = mean(firing_diffs{1}{2}(hinds,:));

figure; hold on; 
for i = 1:length(FDE); 
    
    if (FDE(i,2)<0 && FDL(i,1)>0) 
        
        plot(FDE(i,:)',[1 1]'*FDLA(i),'r');
        plot([1 1]'*FDEA(i),FDL(i,:)','r');
        plot(FDEA(i),FDLA(i),'r.');
        
    else
        
        plot(FDE(i,:)',[1 1]'*FDLA(i),'b');
        plot([1 1]'*FDEA(i),FDL(i,:)','b');
        plot(FDEA(i),FDLA(i),'b.');
    end
end