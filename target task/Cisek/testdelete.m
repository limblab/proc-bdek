switches = [0; find(diff(tt(:,14))~=0); size(tt,1)];
alldays = struct;
for i = 1:(length(switches)-1)
    blockinds = (switches(i)+1):switches(i+1);
    alldays(i).tt = tt(blockinds,:);
    
    % Identify center-out blocks and fill in block_rats variable
    if round(10000*(10*tt(blockinds(1),14) - round(10*tt(blockinds(1),14))))==1 
        alldays(i).ratios = 10000;
        alldays(i).tt(blockinds,14:15) = 10000;
    else
        alldays(i).ratios = tt(blockinds(1),14:15);
    end
end

%%
 figure; hold on; 
 for i = 1:length(ttC)
     t1 = find(bdfM.pos(:,1)>ttC(i,5),1,'first'); 
     t2 = find(bdfM.pos(:,1)<ttC(i,10),1,'last'); 
     
     plot(bdfM.pos(t1:t2,2),bdfM.pos(t1:t2,3)); 
 end