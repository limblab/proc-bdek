%%
figure; hold on;

for i = 1:length(alldays(1).tt)
    
    is = alldays(1).kin.pos(:,1)>alldays(1).tt(i,6) & alldays(1).kin.pos(:,1)<alldays(1).tt(i,7);
    plot(alldays(1).kin.pos(is,2),alldays(1).kin.pos(is,3));
    plot(6.5*cos(alldays(1).tt(i,2)),6.5*sin(alldays(1).tt(i,2)),'.','MarkerSize',5);
end
axis equal; axis square; axis off; 
%%
figure; hold on;

for i = 1:length(alldays(2).tt)
    
    is = alldays(1).kin.pos(:,1)>alldays(2).tt(i,6) & alldays(1).kin.pos(:,1)<alldays(2).tt(i,7);
    plot(alldays(1).kin.pos(is,2),alldays(1).kin.pos(is,3));
    
    plot(6.5*cos(alldays(2).tt(i,2)),6.5*sin(alldays(2).tt(i,2)),'.','MarkerSize',5);
end
axis equal; axis square; axis off