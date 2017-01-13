%% Color map
VM = @(x,k,mu) exp(k.*cos(x-mu))./(2.*pi.*besseli(0,k));

inout = [1.5 2];
positions = 0:0.001:2*pi;
scale = 1;
f = figure('Position',[200 50 1200 950]); hold on;
set(gca,'XTick',[],'YTick',[]);
for count = 1:10
mu = 2*pi*rand;

colr = VM(positions,1,mu);
colr = scale.*(colr/max(colr));

plot((inout(1)-0.01)*cos(0:0.01:2*pi),(inout(1)-0.01)*sin(0:0.01:2*pi),'k');
plot((inout(2)+0.01)*cos(0:0.01:2*pi),(inout(2)+0.01)*sin(0:0.01:2*pi),'k');
for i = 1:length(positions)
    
    plot(inout.*cos(positions(i)),inout.*sin(positions(i)),'Linewidth',2,'Color',[1-colr(i) 1-colr(i) 1]);
    
end
pause;
if strcmp(get(f,'currentkey'),'space')
    close;
    break
end

plot(inout.*cos(mu),inout.*sin(mu),'w','Linewidth',5);

pause; cla;
    
end


%% Wobble

mu = 2*pi*rand;
k = 20;
inout = [1.5 2];

slice = zeros(10000,5);

ds = round(1:10:10000);

for j = 1:5
    for i = 1:length(ds)
        slice(ds(i),j) = mu+vonmisrand(k);
    end
end

sliceint = interp1(ds,slice(ds,1),1:10000);
sliceint2 = interp1(ds,slice(ds,2),1:10000);
sliceint3 = interp1(ds,slice(ds,3),1:10000);
sliceint4 = interp1(ds,slice(ds,4),1:10000);
sliceint5 = interp1(ds,slice(ds,5),1:10000);

sms = smooth(sliceint,20);
sms2 = smooth(sliceint2,20);
sms3 = smooth(sliceint3,20);
sms4 = smooth(sliceint4,20);
sms5 = smooth(sliceint5,20);

f = figure('Position',[200 50 1200 950]); hold on;
set(gca,'XTick',[],'YTick',[]);
for i = 1:length(sms);

    plot(inout(1)*cos(0:0.01:2*pi),inout(1)*sin(0:0.01:2*pi),'k');
    plot(inout(2)*cos(0:0.01:2*pi),inout(2)*sin(0:0.01:2*pi),'k');

    plot(inout.*cos(sms(i)),inout.*sin(sms(i)),'Linewidth',2);
    plot(inout.*cos(sms2(i)),inout.*sin(sms2(i)),'Linewidth',2);
    plot(inout.*cos(sms3(i)),inout.*sin(sms3(i)),'Linewidth',2);
    plot(inout.*cos(sms4(i)),inout.*sin(sms4(i)),'Linewidth',2);
    plot(inout.*cos(sms5(i)),inout.*sin(sms5(i)),'Linewidth',2);

    if strcmp(get(f,'currentkey'),'space')
        
        plot(inout.*cos(mu),inout.*sin(mu),'k','Linewidth',5);
        shg;
    end
    
    
    if strcmp(get(f,'currentkey'),'return')
        close;
        break
    end
    pause(0.001); cla;        
end



%% Add and Replace

mu = 2*pi*rand;
k = 5;
nslices = 10;

slice = zeros(1000,1);
inout = [1.5 2];

for j = 1:1000
    slice(j) = mu+vonmisrand(k);
end


f = figure('Position',[200 50 1200 950]); hold on;
set(gca,'XTick',[],'YTick',[]);
for i = 1:length(slice)-nslices;
    
    plot(inout(1)*cos(0:0.01:2*pi),inout(1)*sin(0:0.01:2*pi),'k');
    plot(inout(2)*cos(0:0.01:2*pi),inout(2)*sin(0:0.01:2*pi),'k');
    
    for j = 1:nslices
        plot(inout.*cos(slice(i+j-1)),inout.*sin(slice(i+j-1,1)),'Linewidth',4);%,'Color',[.8 .8 1]);
    end
    pause(.05); cla;
    if strcmp(get(f,'currentkey'),'return')
        close;
        break
    end
end

%% Imprint and Fade

mu = 2*pi*rand;
k = 5;
ltails = 10;
inout = [1.5 2];

slice = zeros(1000,5);

for q = 1:10
    for j = 1:1000
        slice(j,q) = mu+vonmisrand(k);
    end
end

f = figure('Position',[200 50 1200 950]); hold on;
set(gca,'XTick',[],'YTick',[]);
for i = 1:length(slice)-ltails;
    
    plot(inout(1)*cos(0:0.01:2*pi),inout(1)*sin(0:0.01:2*pi),'k');
    plot(inout(2)*cos(0:0.01:2*pi),inout(2)*sin(0:0.01:2*pi),'k');
    
    for j = 1:ltails
        plot(inout.*cos(slice(i+j-1,1)),inout.*sin(slice(i+j-1,1)),'Linewidth',4,'Color',[(ltails-j)/(ltails+1).*[1 1] 1]);
    end
    
    pause(.05); cla;
    if strcmp(get(f,'currentkey'),'return')
        close;
        break
    end
end

%% Replace All

mu = 2*pi*rand;
k = 25;
inout = [1.5 2];

slice = zeros(1000,5);

for q = 1:5
    for j = 1:1000
        slice(j,q) = mu+vonmisrand(k);
    end
end

f = figure('Position',[200 50 1200 950]); hold on;
set(gca,'XTick',[],'YTick',[]);
for i = 1:length(slice);
    
    plot(inout(1)*cos(0:0.01:2*pi),inout(1)*sin(0:0.01:2*pi),'k');
    plot(inout(2)*cos(0:0.01:2*pi),inout(2)*sin(0:0.01:2*pi),'k');
    
    plot(inout.*cos(slice(i,1)),inout.*sin(slice(i,1)),'Linewidth',2);
    plot(inout.*cos(slice(i+1,2)),inout.*sin(slice(i+1,2)),'Linewidth',2);
    plot(inout.*cos(slice(i+2,3)),inout.*sin(slice(i+2,3)),'Linewidth',2);
    plot(inout.*cos(slice(i+3,4)),inout.*sin(slice(i+3,4)),'Linewidth',2);
    plot(inout.*cos(slice(i+4,5)),inout.*sin(slice(i+4,5)),'Linewidth',2);
    
    pause(.01); cla;
    if strcmp(get(f,'currentkey'),'return')
        close;
        break
    end
end

%% Wobble with Tracer

mu = 2*pi*rand;
k = 20;
nslice = 10;
ntails = 1;
inout = [1.5 2];

slice = zeros(10000,5);

ds = round(1:10:10000);

for j = 1:nslice
    for i = 1:length(ds)
        slice(ds(i),j) = mu+vonmisrand(k);
    end
end

sliceint = cell(nslice,1);
sms = cell(nslice,1);
for i = 1:nslice
    sliceint{i} = interp1(ds,slice(ds,i),1:10000);
    sms{i} = smooth(sliceint{i},20);
end

f = figure('Position',[200 50 1200 950]); hold on;
set(gca,'XTick',[],'YTick',[]);
for i = ntails+1:length(sms{1});
    
    plot(inout(1)*cos(0:0.01:2*pi),inout(1)*sin(0:0.01:2*pi),'k');
    plot(inout(2)*cos(0:0.01:2*pi),inout(2)*sin(0:0.01:2*pi),'k');
    
    for j = 1:ntails
        for q = 1:nslice
            plot(inout.*cos(sms{q}(i-j)),inout.*sin(sms{q}(i-j)),'Linewidth',3,'Color',[j/(ntails+1).*[1 1] 1]);
        end
    end
    
    pause(0.01); cla;
    if strcmp(get(f,'currentkey'),'return')
        close;
        break
    end
end
%% Add and Replace Patches

mu = 2*pi*rand;
k = 25;
patchwidth = pi/20;
lslice = 5;
inout = [1.5 2];

slice = zeros(1000,1);


for j = 1:1000
    slice(j) = mu+vonmisrand(k);
end


f = figure('Position',[200 50 1200 950]); hold on;
set(gca,'XTick',[],'YTick',[]);
for i = 1:length(slice)-lslice;
    
    plot(inout(1)*cos(0:0.01:2*pi),inout(1)*sin(0:0.01:2*pi),'k');
    plot(inout(2)*cos(0:0.01:2*pi),inout(2)*sin(0:0.01:2*pi),'k');
    
    for j = 1:lslice
        p = ang2poly(slice(i+j-1),patchwidth,inout(1),inout(2));
        patch(p(:,1),p(:,2),p(:,3),'EdgeColor','b','EdgeAlpha',0.5,'FaceColor','b','FaceAlpha',0.5)
        %plot([1 2].*cos(slice(i+j-1)),[1 2].*sin(slice(i+j-1,1)),'Linewidth',4);
    end

    pause(.1); cla;
    if strcmp(get(f,'currentkey'),'return')
        close;
        break
    end
end

%% Color Patch
VM = @(x,k,mu) exp(k.*cos(x-mu))./(2.*pi.*besseli(0,k));
mu = 2*pi*rand;
k = 25;
patchwidth = pi/100;
lslice = 200;
inout = [1.5 2];

slice = zeros(1000,1);

range = linspace(0,2*pi,lslice+1); range(end) = [];
probdist = VM(range,k,mu);

for j = 1:10000
    slice(j) = mu+vonmisrand(k);
end


f = figure('Position',[200 50 1200 950]); hold on;
set(gca,'XTick',[],'YTick',[]);
for i = 1:length(slice)-lslice;
    
    plot(inout(1)*cos(0:0.01:2*pi),inout(1)*sin(0:0.01:2*pi),'k');
    plot(inout(2)*cos(0:0.01:2*pi),inout(2)*sin(0:0.01:2*pi),'k');
    
    for j = 1:lslice
        p = ang2poly(range(j),patchwidth,inout(1),inout(2));
        patch(p(:,1),p(:,2),p(:,3),'EdgeColor','b','EdgeAlpha',(probdist(j)/max(probdist))*rand,'FaceColor','b','FaceAlpha',(probdist(j)/max(probdist))*rand)
        %plot([1 2].*cos(slice(i+j-1)),[1 2].*sin(slice(i+j-1,1)),'Linewidth',4);
    end

    pause(.1); cla;
    if strcmp(get(f,'currentkey'),'return')
        close;
        break
    end
end

%% Patch Wobble with Tracer 

mu = 2*pi*rand;
k = 2;
nslice = 20;
ntails = 1;
patchwidth = pi/20;
inout = [1.5 2];

slice = zeros(10000,5);

ds = round(1:10:10000);

for j = 1:nslice
    for i = 1:length(ds)
        slice(ds(i),j) = mu+vonmisrand(k);
    end
end

sliceint = cell(nslice,1);
sms = cell(nslice,1);
for i = 1:nslice
    sliceint{i} = interp1(ds,slice(ds,i),1:10000);
    sms{i} = smooth(sliceint{i},20);
end

f = figure('Position',[200 50 1200 950]); hold on;
set(gca,'XTick',[],'YTick',[]);
for i = ntails+1:length(sms{1});
    
    plot(inout(1)*cos(0:0.01:2*pi),inout(1)*sin(0:0.01:2*pi),'k');
    plot(inout(2)*cos(0:0.01:2*pi),inout(2)*sin(0:0.01:2*pi),'k');
    
    for j = 1:ntails
        for q = 1:nslice
            
            p = ang2poly(sms{q}(i-j),patchwidth,inout(1),inout(2));
            patch(p(:,1),p(:,2),p(:,3),'EdgeColor','b','EdgeAlpha',0.05,'FaceColor','b','FaceAlpha',0.05)
           
        end
    end
    
    pause(0.001); cla;
    if strcmp(get(f,'currentkey'),'return')
        close;
        break
    end
end

%% Imprint and Fade Patches

mu = 2*pi*rand;
k = 5;
ltails = 10;
patchwidth = pi/20;
inout = [1.5 2]; 

slice = zeros(1000,5);

for q = 1:10
    for j = 1:1000
        slice(j,q) = mu+vonmisrand(k);
    end
end

f = figure('Position',[200 50 1200 950]); hold on;
set(gca,'XTick',[],'YTick',[]);
for i = 1:length(slice)-ltails;
    
    plot(inout(1)*cos(0:0.01:2*pi),inout(1)*sin(0:0.01:2*pi),'k');
    plot(inout(2)*cos(0:0.01:2*pi),inout(2)*sin(0:0.01:2*pi),'k');
    
    for j = 1:ltails
        
        p = ang2poly(slice(i+j-1,1),patchwidth,inout(1),inout(2));
        patch(p(:,1),p(:,2),p(:,3),'EdgeColor','b','EdgeAlpha',j/(ltails+1),'FaceColor','b','FaceAlpha',j/(ltails+1))
        
    end
    
    pause(.05); cla;
    if strcmp(get(f,'currentkey'),'return')
        close;
        break
    end
end

%% Num Static

mu = 2*pi*rand;
k = 20;
numslice = 20;

inout = [1.5 2];

slice = zeros(1000,5);


f = figure('Position',[200 50 1200 950]); hold on;
set(gca,'XTick',[],'YTick',[]);


for i = 1:1000

for q = 1:numslice
    slice(q) = mu+vonmisrand(k);
end

plot(inout(1)*cos(0:0.01:2*pi),inout(1)*sin(0:0.01:2*pi),'k');
plot(inout(2)*cos(0:0.01:2*pi),inout(2)*sin(0:0.01:2*pi),'k');

inout_s = 1.5+0.5*rand(numslice,1);
for j = 1:numslice
    plot(inout_s(j).*cos(slice(j)),inout_s(j).*sin(slice(j)),'b.','MarkerSize',30);
end
pause;
for j = 1:numslice
    plot(inout.*cos(slice(j)),inout.*sin(slice(j)),'r','LineWidth',3);
end


pause; plot(inout.*cos(mu),inout.*sin(mu),'Linewidth',10,'Color','g'); pause; cla;
end
