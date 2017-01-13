alpha = 1;
targs = [7  3];
F = figure('Position',[-1596 83 1587 699]); hold on; 
% pL = subplot(1,4,2); 
% pL = subplot('Position',[.035 .05 .2 .8150]); 
pL = subplot('Position',[.035 .05 .35 .8150]); 
hold on;
xlim(xl*.7);
ylim(yl*.7);
zlim(zl*.7);
view(az,el);
axis square;

% plot3([1 1]*.7*xl(1),.7*yl,[1 1]*.7*zl(1),'Color',.5*[1 1 1],'LineWidth',2);
% plot3(.7*xl,[1 1]*.7*yl(2),[1 1]*.7*zl(1),'Color',.5*[1 1 1],'LineWidth',2);
% plot3([1 1]*.7*xl(1),[1 1]*.7*yl(2),.7*zl,'Color',.5*[1 1 1],'LineWidth',2);

plot3([1 1]*.7*xl(1),.7*yl,[1 1]*.7*zl(1),'Color',.5*[1 1 1],'LineWidth',2);
plot3(.7*xl,[1 1]*.7*yl(2),[1 1]*.7*zl(1),'Color',.5*[1 1 1],'LineWidth',2);
plot3([1 1]*.7*xl(1),[1 1]*.7*yl(2),.7*zl,'Color',.5*[1 1 1],'LineWidth',2);
axis off;

t1 = downsample(DR_b2{targs(1)},5); 
t2 = downsample(DR_b2{targs(2)},5);
mean1 = nanmean(DR_b2{targs(1)},1);
mean2 = nanmean(DR_b2{targs(2)},1);
plot3(t1(:,2),t1(:,3),t1(:,4),'.','Color',blue,'MarkerSize',15); %plot3(mean1(:,2),mean1(:,3),mean1(:,4),'b.','MarkerSize',15);
plot3(t2(:,2),t2(:,3),t2(:,4),'.','Color',red,'MarkerSize',15); %plot3(mean2(:,2),mean2(:,3),mean2(:,4),'r.','MarkerSize',15);
% scatter3(t1(:,2),t1(:,3),t1(:,4),15,blue,'MarkerFaceColor',blue,'MarkerFaceAlpha',alpha,'MarkerEdgeColor','none'); %plot3(mean1(:,2),mean1(:,3),mean1(:,4),'b.','MarkerSize',15);
% scatter3(t2(:,2),t2(:,3),t2(:,4),15,red,'MarkerFaceColor',red,'MarkerFaceAlpha',alpha,'MarkerEdgeColor','none'); %plot3(mean2(:,2),mean2(:,3),mean2(:,4),'r.','MarkerSize',15);

all1 = squeeze(nanmean(NTR(targids{1}==targs(1),:,:),1))';
all2 = squeeze(nanmean(NTR(targids{1}==targs(2),:,:),1))';
% set(gcf,'Color','none','inverthardcopy','off');
% 
% plot3(all1(:,2),all1(:,3),all1(:,4),'b','LineWidth',2);
% plot3(all2(:,2),all2(:,3),all2(:,4),'r','LineWidth',2);
%%
% F = figure; hold on; pL = subplot(1,2,1); hold on; 
xlim(xl*.7);
ylim(yl*.7);
zlim(zl*.7);
view(az,el);
axis square;
axis off;

% pR = subplot(1,2,2); hold on; 
% pR = subplot('Position',[0.24 .110 .3347 .8150]); hold on; 
pR = subplot('Position',[0.39 .110 .6150 .8150]); hold on; 
xax =1:126; xax(101:104) = NaN;
plot(pR,xax,.8*ones(length(xax),1),'k');
plot(pR,[1 1],[.8 2],'k'); 
Y = ys;
% Y = Y-repmat(min(Y,[],1),size(Y,1),1);
Y(101:104,:) = nan;
vertlines = [8 49 63 111];
ylim([.8 2.5]);
axis off;
% set(gca,'XTick',[]);
% set(gca,'YTick',[]);
% set(gca,'ZTick',[]);
plot3(pL,[1 1]*.7*xl(1),.7*yl,[1 1]*.7*zl(1),'Color',.5*[1 1 1],'LineWidth',2);
plot3(pL,.7*xl,[1 1]*.7*yl(2),[1 1]*.7*zl(1),'Color',.5*[1 1 1],'LineWidth',2);
plot3(pL,[1 1]*.7*xl(1),[1 1]*.7*yl(2),.7*zl,'Color',.5*[1 1 1],'LineWidth',2);

% set(F,'Color','none','inverthardcopy','off');
% set(pL,'Color','none');

all1(101:104,:) = NaN;
ht = text(pR,-5,1.025,'inverse distance','FontName','Helvetica','FontSize',24);
set(ht,'rotation',90);

ht2 = text(pL,1347.3, 6096.5,-323.3,'PC1','FontName','Helvetica','FontSize',16,'Color',[.5 .5 .5]);
ht3 = text(pL,1344, 6096.2,-331,'PC2','FontName','Helvetica','FontSize',16,'Color',[.5 .5 .5]);
set(ht3,'rotation',-20);
ht4 = text(pL,1337, 6101.5,-323.1,'PC3','FontName','Helvetica','FontSize',16,'Color',[.5 .5 .5]);
set(ht4,'rotation',90);
for i = 1%:size(all1,1); 
    plot3(pL,all1(1:i,2),all1(1:i,3),all1(1:i,4),'k','LineWidth',2);
    cp = plot3(pL,all1(i,2),all1(i,3),all1(i,4),'k.','MarkerSize',15);
    ap1 = plot3(pL,[all1(i,2),mean1(2)],[all1(i,3),mean1(3)],[all1(i,4),mean1(4)],'Color','k','LineWidth',5);
    bp1 = plot3(pL,[all1(i,2),mean2(2)],[all1(i,3),mean2(3)],[all1(i,4),mean2(4)],'Color','k','LineWidth',5);
    ap = plot3(pL,[all1(i,2),mean1(2)],[all1(i,3),mean1(3)],[all1(i,4),mean1(4)],'Color',blue,'LineWidth',3);
    bp = plot3(pL,[all1(i,2),mean2(2)],[all1(i,3),mean2(3)],[all1(i,4),mean2(4)],'Color',red,'LineWidth',3);
    
    ph1 = plot(pR,i,Y(i,1),'.','Color',blue,'MarkerSize',15); 
    ph2 = plot(pR,i,Y(i,2),'.','Color',red,'MarkerSize',15);
    ph3 = plot(pR,1:i,Y(1:i,1),'Color',blue,'LineWidth',3);
    ph4 = plot(pR,1:i,Y(1:i,2),'Color',red,'LineWidth',3);
    if ismember(i,vertlines)
        plot([i i],[.8 2],'-','Color',.5*[1 1 1],'LineWidth',3);
        
        xoff = 6;
        plot(pR,[i-xoff i+xoff i+xoff i-xoff i-xoff],[2.05 2.05 2.3 2.3 2.05],'Color',[.5 .5 .5],'LineWidth',2);
        if ismember(i,vertlines(1:3))
%             plot(i,2.175,'k.','MarkerSize',20);
%             plot(i,2.175,'y.','MarkerSize',12);
            scatter(i,2.175,100,'k','MarkerFaceColor','k','MarkerEdgeColor','none');
            scatter(i,2.175,50,'y','MarkerFaceColor','y','MarkerEdgeColor','none');
            
        end
        if ismember(i,vertlines([1 3 4]))
%             plot(i+xoff/2,2.1125,'.','Color',blue,'MarkerSize',25);
            scatter(i+xoff/2,2.1125,200,'b','MarkerFaceColor',blue,'MarkerEdgeColor','none');
        end
        if i == vertlines(4)
%             plot([i i+xoff/2-.2],[2.175 2.112],'k','LineWidth',5);
%             plot([i i+xoff/2-.2],[2.175 2.112],'y','LineWidth',2);
%             plot(i+xoff/2-.2,2.112,'y.','MarkerSize',12);
            
            plot([i i+xoff/2-.2],[2.175 2.112],'k','LineWidth',4);
            plot([i i+xoff/2-.2],[2.175 2.112],'y','LineWidth',1);
            scatter(i+xoff/2-.2,2.112,50,'y','MarkerFaceColor','y','MarkerEdgeColor','none');
        end 
    end

    pause(0.05); 
%     print(sprintf('C:/Users/limblab/Desktop/pres_figs/TCMC/GIF/Moviepng/t%d',i),'-dpng');
%     if i < 126
%         delete(ap); delete(ap1);
%         delete(bp); delete(bp1);
%         delete(cp);
%         delete(ph1); delete(ph2);
%         delete(ph3); delete(ph4);
%     end
end


%%
wd = 'C:/Users/limblab/Desktop/pres_figs/TCMC/GIF/Moviepng/';
imageNames = dir(fullfile(wd,'*.png'));
imageNames = {imageNames.name}';
imorder = cellfun(@(x) str2double(x((strfind(x,'t')+1):(strfind(x,'.')-1))),imageNames);

sorted_imagenames = cell(size(imageNames));
for i = 1:length(sorted_imagenames)
    sorted_imagenames{i} = imageNames{(imorder==i)};
end

outputVideo = VideoWriter(fullfile(wd,'lowdvid.avi'));
outputVideo.FrameRate = 10;
open(outputVideo); 

for ii = 1:length(imageNames)
    img = imread(fullfile(wd,sorted_imagenames{ii}));
    writeVideo(outputVideo,img)
end
close(outputVideo);