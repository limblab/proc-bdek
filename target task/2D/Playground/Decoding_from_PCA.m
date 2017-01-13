DIR = []; MAG = []; BL = [];
thets = linspace(0,2*pi-pi/4,8)';
for i = 1:length(Dtrial.PMd) % block
    for j = 1:size(Dtrial.PMd{i},1) % trial
        for k = 1:size(Dtrial.PMd{i},2) % epoch
            for m = 1:size(Dtrial.PMd{i}{j,k},2) %time
                y = circshift(Dtrial.PMd{i}{j,k}(:,m),targids{i}(j)-5,1);
                betsPMd = [ones(8,1) cos(thets) sin(thets)]\y;
                DIR.PMd{i,:}{j,k}(:,m) = atan2(betsPMd(3),betsPMd(2));
                MAG.PMd{i,:}{j,k}(:,m) = sqrt(betsPMd(2).^2+betsPMd(3).^2);
                BL.PMd{i,:}{j,k}(:,m) = betsPMd(1); 
                
                y = circshift(Dtrial.M1{i}{j,k}(:,m),targids{i}(j)-5,1);
                betsM1 = [ones(8,1) cos(thets) sin(thets)]\y;
                DIR.M1{i,:}{j,k}(:,m) = atan2(betsM1(3),betsM1(2));
                MAG.M1{i,:}{j,k}(:,m) = sqrt(betsM1(2).^2+betsM1(3).^2);
                BL.M1{i,:}{j,k}(:,m) = betsM1(1); 
            end
        end
    end
end
                
%% Decoded dir v. reach direction, all time
area = 'PMd';
block = 1;
allmag = cell2mat(reshape(vertcat(MAG.(area){:}),1,[]));
conf2size = @(x,sizes) sizes(1)+diff(sizes)*(x-min(allmag))./(max(allmag)-min(allmag));
figure; hold on;
for i = 1:size(DIR.(area){block},1)
    y = mod(cell2mat(DIR.(area){block}(i,:))+2*pi,2*pi);
    m = cell2mat(MAG.(area){block}(i,:));
    sz = conf2size(m,[5 20]); 
    for j = 1:length(y)
        plot(TT{block}(i,align_col),y(j),'b.','MarkerSize',sz(j)); 
    end
end

%% Error in decoded direction over time 
area = 'PMd';
block = 1;
trls = 'all';

if strcmp(trls,'all'); trls = 1:size(DIR.(area){block},1); end
allmag = cell2mat(reshape(vertcat(MAG.(area){:}),1,[]));
conf2size = @(x,sizes) sizes(1)+diff(sizes)*(x-min(allmag))./(max(allmag)-min(allmag));
figure; hold on;
% clrs = {[0 1 0],[0 0 1],[1 0 1],[1 0 0]};
clrs = mat2cell(distinguishable_colors(size(DIR.(area){1},2)),ones(size(DIR.(area){1},2),1),3);
for i = 1:length(trls)
    xs = [0 cumsum(cellfun(@length,DIR.(area){block}(trls(i),:)))];
    for k = 1:size(DIR.(area){block},2)
        x = (xs(k)+1):(xs(k+1));
        y = circ_dist(DIR.(area){block}{trls(i),k},TT{block}(trls(i),align_col));
        m = MAG.(area){block}{trls(i),k};
        
        if ~isnan(sum(m))
            sz = 1-conf2size(m,[0 1]); 
            for j = 1:length(y)-1
                clr = clrs{k}; clr(clr==0) = sz(j); 
    %             plot(x(j),y(j),'.','Color',clrs{k},'MarkerSize',sz(j)); 
                plot([x(j) x(j+1)],[y(j) y(j+1)],'-','Color',clr,'LineWidth',2); 
            end
        end
%         plot(x,y,'-','Color',clrs{k},'LineWidth',0.5);
    end
end
ylim([-pi pi]);

%% Decoder magnitude over time
area = 'M1';
block = 2;
trls = rt{2};

if strcmp(trls,'all'); trls = 1:size(DIR.(area){block},1); end
allmag = cell2mat(reshape(vertcat(MAG.(area){:}),1,[]));
conf2size = @(x,sizes) sizes(1)+diff(sizes)*(x-min(allmag))./(max(allmag)-min(allmag));
figure; hold on;
clrs = {[0 1 0],[0 0 1],[1 0 1],[1 0 0]};
for i = 1:length(trls)
    xs = [0 cumsum(cellfun(@length,DIR.(area){block}(trls(i),:)))];
    for k = 1:size(DIR.(area){block},2)
        x = (xs(k)+1):(xs(k+1));
        y = MAG.(area){block}{trls(i),k};
 
        for j = 1:length(y)-1
%             plot(x(j),y(j),'.','Color',clrs{k},'MarkerSize',sz(j)); 
            plot([x(j) x(j+1)],[y(j) y(j+1)],'-','Color',clrs{k},'LineWidth',1); 
        end
%         plot(x,y,'-','Color',clrs{k},'LineWidth',0.5);
    end
end
% ylim([-pi pi]);

%% Plot reach v decode
block = 2;
area = 'PMd';
regions = {3, 'all'};
dircol = 19;

conf2size = @(x,sizes) sizes(1)+diff(sizes)*(x-min(x(:)))./(max(x(:))-min(x(:)));

[dec,mag] = deal(zeros(size(DIR.(area){block},1),1));
for i = 1:size(DIR.(area){block},1)

    if strcmp(regions{2},'all')
        dec(i) = circ_mean(DIR.(area){block}{i,regions{1}}');
        mag(i) = nanmean(MAG.(area){block}{i,regions{1}});
    else
        decodecell = DIR.(area){block}{i,regions{1}};
        decdir = circ_mean(decodecell(min(regions{2}(1),size(decodecell,2)):min(regions{2}(2),size(decodecell,2)))');
        magcell = MAG.(area){block}{i,regions{1}};
        decmag = nanmean(magcell(min(regions{2}(1),size(magcell,2)):min(regions{2}(2),size(magcell,2))));
        if isempty(decdir)
            dec(i) = NaN;
            mag(i) = NaN;
        else
            dec(i) = decdir;
            mag(i) = decmag;
        end
    end

end
%
dec = mod(dec+2*pi,2*pi);
decsize = conf2size(mag,[5 20]); decsize(isnan(decsize))=0.01;
figure; hold on; 

for i = 1:size(dec,1)
    plot(mod(TT{block}(i,dircol),2*pi),dec(i),'b.','MarkerSize',decsize(i));
end
axis equal
xlim([0 2*pi]); ylim([0 2*pi]); plot([0 2*pi],[0 2*pi],'k');
ylabel('decoded direction','FontSize',14); xlabel('reach direction','FontSize',14);
pause(0.25);  

%% Plot movie of reach v decode
block = 2;
area = 'M1';
% regions = {2, [1  10];...
%            2, [11 20];...
%            2, [21 30];...
%            2, [31 40]};
       
regions = {2, [1   8];...
           2, [9  20];...
           2, [21 30];...
           2, [31 40];...
           3, [1  10];...
           3, [11 20]};
dircol = 10;

conf2size = @(x,sizes) sizes(1)+diff(sizes)*(x-min(x(:)))./(max(x(:))-min(x(:)));

[dec,mag] = deal(zeros(size(DIR.(area){block},1),size(regions,1)));
for i = 1:size(DIR.(area){block},1)
    for j = 1:size(regions,1)
        decodecell = DIR.(area){block}{i,regions{j,1}};
        decdir = circ_mean(decodecell(min(regions{j,2}(1),size(decodecell,2)):min(regions{j,2}(2),size(decodecell,2)))');
        magcell = MAG.(area){block}{i,regions{j,1}};
        decmag = nanmean(magcell(min(regions{j,2}(1),size(magcell,2)):min(regions{j,2}(2),size(magcell,2))));
        if isempty(decdir)
            dec(i,j) = NaN;
            mag(i,j) = NaN;
        else
            dec(i,j) = decdir;
            mag(i,j) = decmag;
        end
    end
end
%
dec = mod(dec+2*pi,2*pi);
decsize = conf2size(mag,[5 20]); decsize(isnan(decsize))=0.01;
figure; hold on; 
for j = 1:size(dec,2)
    cla;
    for i = 1:size(dec,1)
        plot(mod(TT{block}(i,dircol)+2*pi,2*pi),dec(i,j),'b.','MarkerSize',decsize(i,j));
    end
    axis equal
    xlim([0 2*pi]); ylim([0 2*pi]); plot([0 2*pi],[0 2*pi],'k');
    ylabel('decoded direction','FontSize',14); xlabel('reach direction','FontSize',14);
    pause(0.25);  
end

%%
block = 2;
area = 'PMd';
regions = {2, 4;...
           3, 4;...
           4, 4};
dircol = 10;

conf2size = @(x,sizes) sizes(1)+diff(sizes)*(x-min(x(:)))./(max(x(:))-min(x(:)));

[dec,mag] = deal(zeros(size(DIR.(area){block},1),size(regions,1)));
for i = 1:size(DIR.(area){block},1)
    for j = 1:size(regions,1)
        decodecell = DIR.(area){block}{i,regions{j,1}};
        if regions{j,2}(1)>size(decodecell,2) || regions{j,2}(2)
        decdir = circ_mean(decodecell(min(regions{j,2}(1),size(decodecell,2)):min(regions{j,2}(2),size(decodecell,2)))');
        magcell = MAG.(area){block}{i,regions{j,1}};
        decmag = nanmean(magcell(min(regions{j,2}(1),size(magcell,2)):min(regions{j,2}(2),size(magcell,2))));
        if isempty(decdir)
            dec(i,j) = NaN;
            mag(i,j) = NaN;
        else
            dec(i,j) = decdir;
            mag(i,j) = decmag;
        end
    end
end
%
dec = mod(dec+2*pi,2*pi);
decsize = conf2size(mag,[5 20]); decsize(isnan(decsize))=0.01;
figure; hold on; 
for j = 1:size(dec,2)
    for i = 1:size(dec,1)
        plot(TT{block}(i,dircol),dec(i,j),'b.','MarkerSize',decsize(i,j));
    end
    axis equal
    xlim([0 2*pi]); ylim([0 2*pi]); plot([0 2*pi],[0 2*pi],'k');
    ylabel('decoded direction','FontSize',14); xlabel('reach direction','FontSize',14);
    pause(0.25); cla; 
end