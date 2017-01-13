
targs = [7  3];
figure; hold on
xlim(xl*.7);
ylim(yl*.7);
zlim(zl*.7);
view(az,el);
axis square;

plot3([1 1]*.7*xl(1),.7*yl,[1 1]*.7*zl(1),'Color',.5*[1 1 1],'LineWidth',2);
plot3(.7*xl,[1 1]*.7*yl(2),[1 1]*.7*zl(1),'Color',.5*[1 1 1],'LineWidth',2);
plot3([1 1]*.7*xl(1),[1 1]*.7*yl(2),.7*zl,'Color',.5*[1 1 1],'LineWidth',2);
axis off;

t1 = downsample(DR_b2{targs(1)},5); 
t2 = downsample(DR_b2{targs(2)},5);
mean1 = nanmean(DR_b2{targs(1)},1);
mean2 = nanmean(DR_b2{targs(2)},1);
plot3(t1(:,2),t1(:,3),t1(:,4),'.','Color',blue); %plot3(mean1(:,2),mean1(:,3),mean1(:,4),'b.','MarkerSize',15);
plot3(t2(:,2),t2(:,3),t2(:,4),'.','Color',red); %plot3(mean2(:,2),mean2(:,3),mean2(:,4),'r.','MarkerSize',15);

all1 = squeeze(nanmean(NTR(targids{1}==targs(1),:,:),1))';
all2 = squeeze(nanmean(NTR(targids{1}==targs(2),:,:),1))';


% plot3(all1(:,2),all1(:,3),all1(:,4),'b','LineWidth',2);
% plot3(all2(:,2),all2(:,3),all2(:,4),'r','LineWidth',2);
%%
figure; hold on; 
xlim(xl*.7);
ylim(yl*.7);
zlim(zl*.7);
view(az,el);
axis square;
% set(gca,'XTick',[]);
% set(gca,'YTick',[]);
% set(gca,'ZTick',[]);
plot3([1 1]*.7*xl(1),.7*yl,[1 1]*.7*zl(1),'Color',.5*[1 1 1],'LineWidth',2);
plot3(.7*xl,[1 1]*.7*yl(2),[1 1]*.7*zl(1),'Color',.5*[1 1 1],'LineWidth',2);
plot3([1 1]*.7*xl(1),[1 1]*.7*yl(2),.7*zl,'Color',.5*[1 1 1],'LineWidth',2);
axis off;
set(gcf,'Color','none');

for i = 1:size(all1,1); 
    plot3(all1(1:i,2),all1(1:i,3),all1(1:i,4),'k','LineWidth',2);
    cp = plot3(all1(i,2),all1(i,3),all1(i,4),'k.','MarkerSize',15);
    ap1 = plot3([all1(i,2),mean1(2)],[all1(i,3),mean1(3)],[all1(i,4),mean1(4)],'Color','k','LineWidth',5);
    bp1 = plot3([all1(i,2),mean2(2)],[all1(i,3),mean2(3)],[all1(i,4),mean2(4)],'Color','k','LineWidth',5);
    ap = plot3([all1(i,2),mean1(2)],[all1(i,3),mean1(3)],[all1(i,4),mean1(4)],'Color',blue,'LineWidth',3);
    bp = plot3([all1(i,2),mean2(2)],[all1(i,3),mean2(3)],[all1(i,4),mean2(4)],'Color',red,'LineWidth',3);
    pause(0.1); print(sprintf('C:/Users/limblab/Desktop/pres_figs/TCMC/GIF/t%d',i),'-dpng');
    delete(ap); delete(ap1);
    delete(bp); delete(bp1);
    delete(cp);
end