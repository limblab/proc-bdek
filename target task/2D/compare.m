%%
times = (time_bins(2:end)+time_bins(1:end-1))./2;

figure; hold on; 

errorbar(times,nanmean(pc2,2),...
    1.96*nanstd(pc2,[],2)./sqrt(size(pc2,2)),'b');

errorbar(times,nanmean(pc3,2),...
    1.96*nanstd(pc3,[],2)./sqrt(size(pc3,2)),'r');

%%
figure; hold on; 

conc2 = 1 - pv2;
conc3 = 1 - pv3;

errorbar(times,nanmean(conc2,2),...
    1.96*nanstd(conc2,[],2)./sqrt(size(conc2,2)),'b');

errorbar(times,nanmean(conc3,2),...
    1.96*nanstd(conc3,[],2)./sqrt(size(conc3,2)),'r');

%% Mean
bins = 0:0.1:2*pi;
lbin = length(bins);
n = zeros(size(pm2,1),lbin);
x = zeros(size(pm2,1),lbin);
n2 = zeros(size(pm3,1),lbin);
x2 = zeros(size(pm3,1),lbin);
for i = 1:size(pm2,1)
    [n(i,:),x(i,:)] = hist(pm2(i,:),lbin);
    [n2(i,:),x2(i,:)] = hist(pm3(i,:),lbin);
end
%figure; surf((repmat(1:size(pm2,1),36,1))',x,n);
figure; imagesc(bins,fliplr(times),flipud(n));
figure; imagesc(bins,fliplr(times),flipud(n2));
%set(gca,'XTick',bins);
%set(gca,'YTick',times);

%% Var
vbins = 0:0.05:1;
lvbin = length(vbins);
lbin = length(bins);
n = zeros(size(pv2,1),lvbin);
x = zeros(size(pv2,1),lvbin);
n2 = zeros(size(pv3,1),lvbin);
x2 = zeros(size(pv3,1),lvbin);
for i = 1:size(pv2,1)
    [n(i,:),x(i,:)] = hist(pv2(i,:),lvbin);
    [n2(i,:),x2(i,:)] = hist(pv3(i,:),lvbin);
end
%figure; surf((repmat(1:size(pm2,1),36,1))',x,n);
figure; imagesc(vbins,fliplr(times),flipud(n));
figure; imagesc(vbins,fliplr(times),flipud(n2));
%set(gca,'XTick',bins);
%set(gca,'YTick',times);

%%
av_code = zeros(length(C2),size(C2{1},2));
av_code_cos = zeros(length(C2C),size(C2C{1},2));
for i = 1:length(C2) 
    av_code(i,:) = nanmean(C2{i},1);
    av_code_cos(i,:) = nanmean(C2C{i},1);
end
%figure; imagesc(wrapped_cents,times,av_code);
figure; imagesc(wrapped_cents,times,av_code_cos);

%%
av_code = zeros(length(C3),size(C3{1},2));
av_code_cos2 = zeros(length(C3C),size(C3C{1},2));
for i = 1:length(C3) 
    av_code(i,:) = nanmean(C3{i},1);
    av_code_cos2(i,:) = nanmean(C3C{i},1);
end
%figure; imagesc(wrapped_cents,times,av_code);
figure; imagesc(wrapped_cents,times,av_code_cos2);


%% Booted errorbar.
m2 = mean(booted2,1);
m3 = mean(booted3,1);

e_mat2 = sort(booted2);
e_mat3 = sort(booted3);
l_e_ind = round(size(booted2,1)*0.025);
h_e_ind = size(booted2,1) - l_e_ind + 1;

figure; hold on;
patch([times fliplr(times)],[e_mat2(l_e_ind,:) fliplr(e_mat2(h_e_ind,:))],'b','FaceAlpha',0.5,'EdgeAlpha',0.5);
plot(times,m2,'b');

patch([times fliplr(times)],[e_mat3(l_e_ind,:) fliplr(e_mat3(h_e_ind,:))],'r','FaceAlpha',0.5,'EdgeAlpha',0.5);
plot(times,m3,'r');

