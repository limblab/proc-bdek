function [kinstruct] = kinbytrial(alldays)

TS = alldays(1).kin.pos(:,1);
kinstruct = cell(length(alldays),1);
for i = 1:length(alldays)
    
    for j = 1:size(alldays(i).tt,1)
        
        clc;fprintf('%d/%d - %d/%d\n',i,length(alldays),j,size(alldays(i).tt,1));
        
        t1 = alldays(i).tt(j,5)-0.25;
        t2 = alldays(i).tt(j,5);
        t3 = alldays(i).tt(j,6);
        t4 = alldays(i).tt(j,7);
        kinstruct{i}{j,:}.pos{1,:} = alldays(1).kin.pos(TS>t1&TS<t2,2:3);
        kinstruct{i}{j,:}.pos{2,:} = alldays(1).kin.pos(TS>t2&TS<t3,2:3);
        kinstruct{i}{j,:}.pos{3,:} = alldays(1).kin.pos(TS>t3&TS<t4,2:3);
       
    end
end
        
        
        