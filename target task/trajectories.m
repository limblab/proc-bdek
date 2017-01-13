tt = alldays(1).tt;
waits = tt(:,[6 7]);
figure; hold on;
[spos,epos] = deal(zeros(length(waits),2));
for i = 1:length(waits)
    startindw = find(BDF.pos(:,1)>waits(i,1),1,'first');
    endindw = find(BDF.pos(:,1)<waits(i,2),1,'last');
    x = BDF.pos(startindw:endindw,2);
    y = BDF.pos(startindw:endindw,3);
%     plot(x,y,'g','LineWidth',5);
    r = sqrt(x(end).^2+y(end).^2);
    
    spos(i,:) = [x(1) y(1)];
    epos(i,:) = [x(end) y(end)];
%     plot(r.*cos(0:0.01:2*pi),r.*sin(0:0.01:2*pi),'k');
%     plot(r.*cos(alldays(2).tt(i,10)),r.*sin(alldays(2).tt(i,10)),'r.','MarkerSize',14);
    %pause;
    %cla;

end
%%
figure; hold on;
for i = 1:length(moves)
    startindw = find(BDF.pos(:,1)>waits(i,1),1,'first');
    endindw = find(BDF.pos(:,1)<waits(i,2),1,'last');
    
    plot(BDF.vel(startindw:endindw,1)-BDF.vel(startindw,1),BDF.vel(startindw:endindw,3),'g');
    
    startind = find(BDF.pos(:,1)>moves(i,1),1,'first');
    endind = find(BDF.pos(:,1)<moves(i,2),1,'last');
    
    plot(BDF.vel(startind:endind,1)-BDF.vel(startindw,1),BDF.vel(startind:endind,3),'b');
    
    drawnow; pause;

end