%%
com1 = 1; 
com2 = 2;
[LIK,LIKt,LIKs,LIKc,pcs,scs] = deal(cell(length(Rec.all{1}),1));
for i = 1:length(Rec.all{1})
  
    L = cellfun(@(x) x{i}, Rec.all,'UniformOutput',0);
    Lm = cellfun(@mean,L,'UniformOutput',0);
    LIKt{i} = vertcat(Lm{:});
  
    LIK{i} = LIKt{i};
    %%% for use with PPC_Sanger output
%     LIK{i} = LIKt{i}(:,1:(end-1));
%     LIKs{i} = LIKt{i}(:,end);
    %%%
   
    [pcs{i}, scs{i}] = princomp(LIK{i});
end

[err1,err2] = deal(zeros(length(Rec.all),1));
for i = 1:length(Rec.all)
    
    max1 = find(LIK{com1}(i,:)==max(LIK{com1}(i,:)));
    likmax1 = [LIK{com1}(i,max1:end) LIK{com1}(i,1:(max1-1))];
    likalign1 = [likmax1((end/2):end), likmax1(1:(end/2-1))];
    
    err1(i) = ((max1-314)./100)/pi*180;
    
    max2 = find(LIK{com2}(i,:)==max(LIK{com2}(i,:)));
    likmax2 = [LIK{com2}(i,max2:end) LIK{com2}(i,1:(max2-1))];
    likalign2 = [likmax2((end/2):end), likmax2(1:(end/2-1))];
    
    err2(i) = ((max2-314)./100)/pi*180;
    
    LIKc{com1}(i,:) = likalign1;
    LIKc{com2}(i,:) = likalign2;
    
end

%LIK = LIK;

%%    
minv = min(cellfun(@(x) min(x(:)),LIKt));
maxv = max(cellfun(@(x) max(x(:)),LIKt));

normf = @(x,mn,mx) 64*(x-mn)./(mx-mn); 

midtimes = 0.5*(time_bins(2:end)+time_bins(1:end-1));
yinterp = zeros(1000,628);
figure; 
for i = 1:length(LIK)
    subplot(length(LIK)+1,1,i); hold on;
%     surf(repmat((midtimes)',1,628),repmat(linspace(-180,180,628),size(LIK{i},1),1),...
%         normf(LIK{i},minv,maxv),normf(LIK{i},minv,maxv),'EdgeColor','none','CDataMapping','direct','FaceColor','interp');
%     view(0,90); 
    
    for j = 1:size(LIK{1},2)
        yinterp(:,j) = interp1(midtimes',normf(LIK{i}(:,j),minv,maxv),linspace(midtimes(1),midtimes(end),1000)');
    end
    xs = linspace(-180,180,628);
    ys = linspace(midtimes(1),midtimes(end),1000);
    
    image(ys,xs,yinterp');
    
%     surf([midtimes' midtimes'],[-195*ones(10,1) -180*ones(10,1)],...
%         repmat(normf(1-LIKs{i},minv,maxv),1,2),repmat(normf(1-LIKs{i},minv,maxv),1,2),...
%         'EdgeColor','none','CDataMapping','direct','FaceColor','interp');
    
   % set(gcf,'Renderer','Zbuffer');
    
    xlim([midtimes(1) midtimes(end)]); if i ==(length(LIK)+1); xlabel('Time from Target On (ms)','FontSize',14); end
    ylim([-180 180]); if i ==2; ylabel('Distance from Reach Direction','FontSize',14); end
    %shading interp;
end

diff_rec = (LIK{com2}-LIK{com1});%matchrec;
diff_rec_n = normf(diff_rec,min(diff_rec(:)),max(diff_rec(:)));

subplot(length(LIK)+1,1,length(LIK)+1); 
for j = 1:size(LIK{1},2)
    yinterp(:,j) = interp1(midtimes',diff_rec_n(:,j),linspace(midtimes(1),midtimes(end),1000)');
end
xs = linspace(-180,180,628);
ys = linspace(midtimes(1),midtimes(end),1000);
    
image(ys,xs,yinterp');


%% Bootstrap
[contain0,pcontain0] = deal(zeros(length(Rec.all),628));
for i = 1:length(Rec.all)
    
    Rt = Rec.all{i};
    
    [~,~,~,R1rand] = boot_bounds(1000,@mean,Rt{com1},2.5,97.5);
    [~,~,~,R2rand] = boot_bounds(1000,@mean,Rt{com2},2.5,97.5);

    base_diff = repmat(mean(R2rand,2) - mean(R1rand,2),1,628);
    ds = R2rand - (R1rand + base_diff);
    
    for j = 1:size(ds,2)
        
        ranked_d = sortrows(ds(:,j));
        contain0(i,j) = ranked_d(round(0.025*1000))*ranked_d(round(0.975*1000)) < 0;
        
        pval = -(find(ranked_d>0,1,'first') - 500)/500;
        
        if isempty(pval); 
            pcontain0(i,j) = sign(ranked_d(1)); 
        else
            pcontain0(i,j) = pval;
        end
        
    end
end

sigposneg = pcontain0;
sigposneg(pcontain0 > 0.975)=1; sigposneg(pcontain0 < - 0.975)= -1; 
sigposneg(~ismember(sigposneg,[-1 1])) = 0;

% spninterp = zeros(1000,628);
% for j = 1:size(LIK{1},2)
%     spninterp(:,j) = interp1(midtimes',sigposneg(:,j),linspace(midtimes(1),midtimes(end),1000)');
% end
%%
figure; hold on; subplot(2,1,1); hold on; 
imagesc(ys,xs,sigposneg');%,'FaceColor','interp'); 
xlim([midtimes(1) midtimes(end)]); ylim([-180 180]); 
view(0,90); 
%%
minmax = @(x) [min(x) max(x)];
vmfunc = @(a) a(1) + a(2)*exp(a(3)*cos(wrapped_cents-a(4)));
vmfunc2 = @(a) a(1) + a(2)*exp(a(3)*cos(wrapped_cents-a(4))) + a(5)*exp(a(6)*-cos(wrapped_cents-a(4)));
fwhm = @(x) sum((x - 0.5*(min(x) + max(x))) > 0)*180/(pi*100);


figure; hold on; 

[a,k,wid] = deal(zeros(1,size(LIK{com1},1)));
[matchrec] = deal(zeros(size(LIK{com1})));
for i = 1:size(LIK{com1},1)
    
    p = polyfit(minmax(LIK{com1}(i,:)),minmax(LIK{com2}(i,:)),1);
    %p = polyfit(LIK{com1}(i,:),LIK{com2}(i,:),1);  
    a(i) = p(1); k(i) = p(2);
    
    wid(i) = fwhm(LIK{com2}(i,:)) - fwhm(LIK{com1}(i,:));
%     A1(i,:) = fminsearch(@(x) sum((vmfunc(x)-LIK{com1}(i,:)).^2),[0 1 1 pi]);
%     A2(i,:) = fminsearch(@(x) sum((vmfunc(x)-LIK{com2}(i,:)).^2),[0 1 1 pi]);
%     A1(i,:) = fminsearch(@(x) sum((vmfunc2(x)-LIK{com1}(i,:)).^2),[0 1 1 pi 0 1]);
%     A2(i,:) = fminsearch(@(x) sum((vmfunc2(x)-LIK{com2}(i,:)).^2),[0 1 1 pi 0 1]);
% 
% %     
%     fit1(i,:) = vmfunc2(A1(i,:));
%     fit2(i,:) = vmfunc2(A2(i,:));
%     
%     k(i) = diff(minmax(LIK{com2}(i,:)))/diff(minmax(LIK{com1}(i,:)));
%     a(i) = median(LIK{com2}(i,:))-median(LIK{com1}(i,:));
    
    %base_diff = mean(LIK{com2},2) - mean(LIK{com1},2);
%     base_diff = max(LIK{com2},[],2) - max(LIK{com1},[],2);
    
    maxdiff = max(LIK{com2},[],2) - max(LIK{com1},[],2);
    mindiff = min(LIK{com2},[],2) - min(LIK{com1},[],2);
%     meddiff = median(LIK{com2},2) - median(LIK{com1},2);
%     
%     gain_diff(i) = diff(minmax(LIK{com2}(i,:)))./diff(minmax(LIK{com1}(i,:)));
    
%     base_diff_rep = repmat(base_diff,1,628);
%     gain_diff_rep(i,:) = repmat(gain_diff(i),1,628);
%     
%    matchrec = base_diff_rep + LIK{com1};
    
    matchrec(i,:) = polyval(p,LIK{com1}(i,:));
%     matchrec(i,:) = k(i)*(LIK{com1}(i,:)-min(LIK{com1}(i,:))) + a(i);
%     diff_rec(i,:) = LIK{com2}(i,:) - polyval(p,LIK{com1}(i,:));
    
    
end

%matchrec = (LIK{com1}-repmat(min(LIK{com1},[],2),1,628)).*gain_diff_rep + repmat(min(LIK{com2},[],2),1,628);

diff_rec = (LIK{com2}-LIK{com1});%matchrec;

%diff_rec_n = normf(diff_rec,min(LIKt{2}(:)),max(LIKt{2}(:))) + 32-median(median(normf(diff_rec,min(LIKt{2}(:)),max(LIKt{2}(:)))));
diff_rec_n = normf(diff_rec,min(diff_rec(:)),max(diff_rec(:)));

subplot(4,1,1); 

for j = 1:size(LIK{1},2)
    yinterp(:,j) = interp1(midtimes',diff_rec_n(:,j),linspace(midtimes(1),midtimes(end),1000)');
end
xs = linspace(-180,180,628);
ys = linspace(midtimes(1),midtimes(end),1000);
    
image(ys,xs,yinterp');

% surf(repmat((midtimes)',1,628),repmat(linspace(-180,180,628),size(LIK{1},1),1),...
%     ones(size(LIK{1})), diff_rec_n,'EdgeColor','none','FaceColor','interp'); 
% 
% xlim([midtimes(1) midtimes(end)]); xlabel('Time from Target On (ms)','FontSize',14);
% ylim([-180 180]); %ylabel('Distance from Reach Direction','FontSize',14);
% view(0,90); 
%colorbar;

% subplot(2,1,2);
% plot(midtimes,base_diff); ylabel('gain','FontSize',14);

subplot(4,1,2);
plot(midtimes,mindiff); ylabel('dmin','FontSize',14);

subplot(4,1,3);
plot(midtimes,maxdiff); ylabel('dmax','FontSize',14);

subplot(4,1,4);
plot(midtimes,wid); ylabel('width','FontSize',14);

%% Average-matched plots
figure; hold on;

subplot(2,1,1); hold on;
surf(repmat((midtimes)',1,628),repmat(linspace(-180,180,628),size(LIK{com1},1),1),...
    matchrec,matchrec,'EdgeColor','none','FaceColor','interp');

xlim([midtimes(1) midtimes(end)]); xlabel('Time from Target On (ms)','FontSize',14);
ylim([-180 180]); %ylabel('Distance from Reach Direction','FontSize',14);
view(0,90); 

subplot(2,1,2); hold on; 
surf(repmat((midtimes)',1,628),repmat(linspace(-180,180,628),size(LIK{com1},1),1),...
    normf(LIK{com2},minv,maxv),normf(LIK{com2},minv,maxv),'EdgeColor','none','CDataMapping','direct','FaceColor','interp');

xlim([midtimes(1) midtimes(end)]); xlabel('Time from Target On (ms)','FontSize',14);
ylim([-180 180]); %ylabel('Distance from Reach Direction','FontSize',14);
view(0,90); 

set(gcf,'Renderer','Zbuffer');
    