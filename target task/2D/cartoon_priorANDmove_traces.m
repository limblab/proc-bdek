binwidth = pi/50;

for i = 1:10000
    
    dlow(i) = pi/2 + vonmisrand(50);
    dhi(i) = pi/2 + vonmisrand(5);
  
end

    inneredge = 1000;
    
    [nl,xl] = hist(dlow,0:binwidth:2*pi);
    [nh,xh] = hist(dhi,0:binwidth:2*pi);
   
figure; hold on;
for i = 1:length(xl)
    
    p = ang2poly(xl(i),binwidth+.007169,inneredge,inneredge + nl(i));
    
%     plx(:,i) = p(:,1); 
%     ply(:,i) = p(:,2);
    
    patch(p(:,1),p(:,2),'b');
end
axis equal
  
figure; hold on ;
for i = 1:length(xh)
    p = ang2poly(xh(i),binwidth+.007169,inneredge,inneredge + nh(i));
    
%     phx(:,i) = p(:,1); 
%     phy(:,i) = p(:,2);
%     
    
    patch(p(:,1),p(:,2),'r');
end

axis equal

%%
day = 2;
tracex = nan(length(alldays(day).tt),1000*ceil(alldays(day).tt(:,7)-alldays(day).tt(:,6)));
tracey = tracex;
for i = 1:length(alldays(day).tt)
    
    movestart = find(alldays(day).kin.pos(:,1) > alldays(day).tt(i,6),1,'first');
    moveend = find(alldays(day).kin.pos(:,1) < alldays(day).tt(i,7),1,'last');
    
    move = alldays(day).kin.pos(movestart:moveend,2:3)';
   
    tracex(i,1:length(move)) = move(1,:);
    tracey(i,1:length(move)) = move(2,:);
    
    clc
    fprintf('%d/%d\n',i,length(alldays(day).tt));
    
end
figure; plot(tracex' - 1,tracey' + 0.5,'b');
p = ang2poly(0,2*pi,8,10);
hold on; 
patch(p(:,1),p(:,2),'w');

