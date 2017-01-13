figure; hold on; 

subplot(2,4,1);
imagesc(nanmean(M.pmd{1}.d(:,:,linds{1}),3));
subplot(2,4,5);
imagesc(nanmean(M.pmd{1}.d(:,:,hinds{1}),3));
cc = caxis;
for i = 1:3
    
    subplot(2,4,i+1);
    imagesc(nanmean(M.pmd{i+1}.d(:,:,linds{i+1}),3),cc);
    subplot(2,4,i+5);
    imagesc(nanmean(M.pmd{i+1}.d(:,:,hinds{i+1}),3),cc);
    
end

%%
figure; hold on; 

subplot(2,4,1);
imagesc(nanmean(M.m1{1}.d(:,:,linds{1}),3));
subplot(2,4,5);
imagesc(nanmean(M.m1{1}.d(:,:,hinds{1}),3));
cc = caxis;
for i = 1:3
    
    subplot(2,4,i+1);
    imagesc(nanmean(M.m1{i+1}.d(:,:,linds{i+1}),3),cc);
    subplot(2,4,i+5);
    imagesc(nanmean(M.m1{i+1}.d(:,:,hinds{i+1}),3),cc);
    
end

%%
ABzero = cellfun(@(x) x - repmat(nanmean(x(:,1:7,:),2),1,46,1),AB,'Uni',0);
ABmzero = cellfun(@(x) x - repmat(nanmean(x(:,1:7,:),2),1,46,1),ABm,'Uni',0);

%%
figure; hold on; 
for i = 1:length(alldays(2).tt)
    
    p1 = subplot(2,1,1); hold on;
    act = circshift(M.pmd{2}.d(:,:,i),-targids{2}(i),1);
    imagesc([1 size(act,2)],[1 8],act);
    rch = mod(mod(alldays(2).tt(i,10)+4*pi,2*pi)./(pi/4)+1,9);
    slcs = mod(mod(alldays(2).slices(i,:)+4*pi,2*pi)./(pi/4)+1,9);

    plot(repmat([1;size(M.pmd{2}.d,2)],1,5),repmat(slcs,2,1),'r');
    plot([1 size(M.pmd{2}.d,2)],[rch rch],'w','LineWidth',3);
    xlim([1 size(act,2)]); ylim([1 8]);
    
    p2 = subplot(2,1,2); hold on;
    act = circshift(M.m1{2}.d(:,:,i),-targids{2}(i),1);
    imagesc([1 size(act,2)],[1 8],act);

    plot(repmat([1;size(M.m1{2}.d,2)],1,5),repmat(slcs,2,1),'r');
    plot([1 size(M.m1{2}.d,2)],[rch rch],'w','LineWidth',3);
    xlim([1 size(act,2)]); ylim([1 8]);
    pause; cla(p1); cla(p2);
end
%%
figure; hold on; 
for i = 1:length(alldays(2).tt)
    
%     actpmd = diff(M.pmd{2}.d([1 5],:,i)./repmat(nanmean(M.pmd{2}.d([3 7],:,i),1),2,1),[],1);
%     actm1 = diff(M.m1{2}.d([1 5],:,i)./repmat(nanmean(M.m1{2}.d([3 7],:,i),1),2,1),[],1);
   
    actpmd = M.pmd{2}.d(5,:,i)./nanmean(M.pmd{2}.d([3 7],:,i),1);
    actm1 = M.m1{2}.d(5,:,i)./nanmean(M.m1{2}.d([3 7],:,i),1);
    if alldays(2).tt(i,3)>25
         plot(actpmd,'b','LineWidth',3); plot(actm1,'b:','LineWidth',3); 
    else
         plot(actpmd,'r','LineWidth',3); plot(actm1,'r:','LineWidth',3); 
    end
    pause; cla;
end
%%
blck = 4;
figure; hold on; 
subplot(1,2,1); hold on; 
set1 = linds{blck};
plot(nanmean(M.pmd{blck}.d(5,:,set1),3)','b');
plot(nanmean(nanmean(M.pmd{blck}.d([4 6],:,set1),1),3)','m');
plot(nanmean(nanmean(M.pmd{blck}.d([3 7],:,set1),1),3)','r'); 
plot(nanmean(nanmean(M.pmd{blck}.d([2 8],:,set1),1),3)','g');
plot(nanmean(nanmean(M.pmd{blck}.d(1,:,set1),1),3)','k');
yl1 = ylim;

subplot(1,2,2); hold on; 
set2 = hinds{blck};
plot(nanmean(M.pmd{blck}.d(5,:,set2),3)','b');
plot(nanmean(nanmean(M.pmd{blck}.d([4 6],:,set2),1),3)','m');
plot(nanmean(nanmean(M.pmd{blck}.d([3 7],:,set2),1),3)','r'); 
plot(nanmean(nanmean(M.pmd{blck}.d([2 8],:,set2),1),3)','g');
plot(nanmean(nanmean(M.pmd{blck}.d(1,:,set2),1),3)','k');
ylim(yl1);

%%
figure; hold on; 
plot(nanmean(M.m1{2}.d(5,:,linds{2}),3)','b');
plot(nanmean(nanmean(M.m1{2}.d([4 6],:,linds{2}),1),3)','g');
plot(nanmean(nanmean(M.m1{2}.d([3 7],:,linds{2}),1),3)','r'); 
plot(nanmean(nanmean(M.m1{2}.d([2 8],:,linds{2}),1),3)','k');

figure; hold on; 
plot(nanmean(M.m1{2}.d(5,:,hinds{2}),3)','b');
plot(nanmean(nanmean(M.m1{2}.d([4 6],:,hinds{2}),1),3)','g');
plot(nanmean(nanmean(M.m1{2}.d([3 7],:,hinds{2}),1),3)','r'); 
plot(nanmean(nanmean(M.m1{2}.d([2 8],:,hinds{2}),1),3)','k');

%%
Mshift = [];
for blck = 1:length(M.pmd)
    for i = 1:size(M.pmd{blck}.d,3)
        Mshift.pmd{blck}.d(:,:,i) = circshift(M.pmd{blck}.d(:,:,i),targids{blck}(i)-5,1); 
        Mshift.m1{blck}.d(:,:,i) = circshift(M.m1{blck}.d(:,:,i),targids{blck}(i)-5,1);
    end
end

%% 
block = 2;
area = 'pmd';
region = 25:35;

[decode,conf,maxrep,minrep,offsetlev] = deal(zeros(size(Mshift.(area){block}.d,3),size(Mshift.(area){block}.d,2)));
for i = 1:size(Mshift.(area){block}.d,3)
    for j = 1:size(Mshift.(area){block}.d,2)
        xs = (0:pi/4:(2*pi-pi/4))';
        ys = Mshift.(area){block}.d(:,j,i);
        
        bs = [ones(8,1) cos(xs) sin(xs)]\ys;
        decode(i,j) = atan2(bs(3),bs(2));
        conf(i,j) = sqrt(bs(2).^2+bs(3).^2);
        
        maxrep(i,j) = bs(1)+.5*conf(i,j);
        minrep(i,j) = bs(1)-0.5*conf(i,j);
        offsetlev(i,j) = bs(1);
    end
end
conf2size = @(x,sizes) sizes(1)+diff(sizes)*(x-min(x(:)))./(max(x(:))-min(x(:)));
figure; hold on; 
decdir = mod(circ_mean(decode(:,region),[],2)+2*pi,2*pi);
decconf = nanmean(conf(:,region),2);
decsize = conf2size(decconf,[5 20]); decsize(isnan(decsize))=1;
for i = 1:length(linds{block})
    plot(alldays(block).tt(linds{block}(i),10),decdir(linds{block}(i)),'b.','MarkerSize',decsize(linds{block}(i)));
end
for i = 1:length(hinds{block})
    plot(alldays(block).tt(hinds{block}(i),10),decdir(hinds{block}(i)),'r.','MarkerSize',decsize(hinds{block}(i)));
end
axis equal
xlim([0 2*pi]); ylim([0 2*pi]); plot([0 2*pi],[0 2*pi],'k');
        
%%
block = 2;
area = 'pmd';
region = 25:35;
dircol = 10;

[decode,conf,maxrep,minrep,offsetlev] = deal(zeros(size(Mshift.(area){block}.d,3),1));
[decodeall,confall,maxrepall,minrepall,offsetlevall] = ...
    deal(zeros(size(Mshift.(area){block}.d,3),size(Mshift.(area){block}.d,2)));
Bcos = zeros(size(Mshift.(area){block}.d,3),100);
for i = 1:size(Mshift.(area){block}.d,3)
    xs = repmat((0:pi/4:(2*pi-pi/4))',1,length(region));
    ys = Mshift.(area){block}.d(:,region,i);
    xslong = linspace(0,2*pi,100)';
    
    bs = [ones(numel(xs),1) cos(reshape(xs,[],1)) sin(reshape(xs,[],1))]\reshape(ys,[],1);
    bcosine = [ones(numel(xslong),1) cos(xslong) sin(xslong)]*bs;
    Bcos(i,:) = circshift(bcosine,50-round(alldays(block).tt(i,dircol)*50/pi));
    
%     bs = glmfit([cos(reshape(xs,[],1)) sin(reshape(xs,[],1))],reshape(ys,[],1),'poisson');
%     bcosine = glmval(bs,[cos(xslong) sin(xslong)],'log');
%     Bcos(i,:) = circshift(bcosine,50-round(alldays(block).tt(i,dircol)*50/pi));
    
    decode(i,:) = atan2(bs(3),bs(2));
    conf(i,:) = sqrt(bs(2).^2+bs(3).^2);
    maxrep(i,:) = bs(1)+.5*conf(i,:);
    minrep(i,:) = bs(1)-.5*conf(i,:);
    offsetlev(i,:) = bs(1);
    
    for j = 1:size(Mshift.(area){block}.d,2)
        xs = (0:pi/4:(2*pi-pi/4))';
        ys = Mshift.(area){block}.d(:,j,i);
        
        bs = [ones(8,1) cos(xs) sin(xs)]\ys;
        decodeall(i,j) = atan2(bs(3),bs(2));
        confall(i,j) = sqrt(bs(2).^2+bs(3).^2);
        
        maxrepall(i,j) = bs(1)+.5*confall(i,j);
        minrepall(i,j) = bs(1)-0.5*confall(i,j);
        offsetlevall(i,j) = bs(1);
    end
end
conf2size = @(x,sizes) sizes(1)+diff(sizes)*(x-min(x(:)))./(max(x(:))-min(x(:)));
decodeall = mod(decodeall+2*pi,2*pi);

figure; hold on; 
decdir = mod(decode+2*pi,2*pi);
decconf = nanmean(conf,2);
decsize = conf2size(decconf,[5 20]); decsize(isnan(decsize))=0.01;
for i = 1:length(linds{block})
    plot(alldays(block).tt(linds{block}(i),dircol),decdir(linds{block}(i)),'b.','MarkerSize',decsize(linds{block}(i)));
end
for i = 1:length(hinds{block})
    plot(alldays(block).tt(hinds{block}(i),dircol),decdir(hinds{block}(i)),'r.','MarkerSize',decsize(hinds{block}(i)));
end
axis equal
xlim([0 2*pi]); ylim([0 2*pi]); plot([0 2*pi],[0 2*pi],'k');
ylabel('decoded direction','FontSize',14); xlabel('reach direction','FontSize',14);
        
%%
sz = conf2size(confall,[5 20]); sz(isnan(sz)) = 1;
figure; hold on; 
for j = 1:size(sz,2)
    cla;
   
    plot([0 2*pi],[0 2*pi],'k');
    for i = 1:length(linds{block})
        plot(alldays(block).tt(linds{block}(i),dircol),mod(decodeall(linds{block}(i),j)+2*pi,2*pi),'b.','MarkerSize',sz(linds{block}(i),j));
    end
    for i = 1:length(hinds{block})
        plot(alldays(block).tt(hinds{block}(i),dircol),mod(decodeall(hinds{block}(i),j)+2*pi,2*pi),'r.','MarkerSize',sz(hinds{block}(i),j));
    end
    axis equal;
    xlim([0 2*pi]); ylim([0 2*pi]); 
    title(sprintf('%d',j));
    pause(0.025); 
end
    
%%
sz = conf2size(confall,[5 20]); sz(isnan(sz)) = 1;
for i = 1:size(sz,1)
    szterp(i,:) = interp1(1:46,sz(i,:),1:0.1:46);
    decodeterp(i,:) = interp1(1:46,decodeall(i,:),1:0.1:46);
end
figure; hold on; 
for j = 1:size(szterp,2)
    cla;
   
    plot([0 2*pi],[0 0],'k');
    plot(alldays(block).tt(linds{block},dircol),circ_dist(decodeterp(linds{block},j),alldays(block).tt(linds{block},dircol)),'b.','MarkerSize',5);
    plot(alldays(block).tt(hinds{block},dircol),circ_dist(decodeterp(hinds{block},j),alldays(block).tt(hinds{block},dircol)),'r.','MarkerSize',5);

    axis equal;
    xlim([0 2*pi]); ylim([-pi pi]); 
    title(sprintf('%d',j));
    pause(0.001); 
end