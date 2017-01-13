AvB1 = {cell2mat(cellfun(@(x) x(:,1),TAB(:,1)','Uni',0)), cell2mat(cellfun(@(x) x(:,2),TAB(:,1)','Uni',0))};
AvB2 = {cell2mat(cellfun(@(x) x(:,1),TAB(:,2)','Uni',0)), cell2mat(cellfun(@(x) x(:,2),TAB(:,2)','Uni',0))};
AvB = {cell2mat([AvB1(:,1) AvB2(:,1)]) , cell2mat([AvB1(:,2) AvB2(:,2)])};
%%
pnAB = zeros(size(AvB{1},1),2);
for i = 1:size(AvB{1},1)
    pnAB(i,1) = ttest(AvB{1}(i,:))*sign(nanmean(AvB{1}(i,:)));
    pnAB(i,2) = ttest(AvB{2}(i,:))*sign(nanmean(AvB{2}(i,:)));
end

simul = find(sum(pnAB,2)==2);
norep = find(sum(abs(pnAB),2)==0);
favorstay = find(pnAB(:,1)==1 & pnAB(:,2)==0);
favorswap = find(pnAB(:,1)==0 & pnAB(:,2)==1);
strongbiasA = find(pnAB(:,1)==1 & pnAB(:,2)==-1);
strongbiasB = find(pnAB(:,1)==-1 & pnAB(:,2)==1);



%%
spans0 = @(x) (x(:,1).*x(:,2) > 0).*sign(mean(x,2));
region = mem;
barea = 'm1';
[favorA,favorB,bothAB,neitherAB,AoverB,AorB,AtoB,BtoA,dbias_time,dbias_time_m1,AoverB_m] = deal(zeros(size(M.(lower(barea)){2}.d,3),1));
for i = 1:size(M.(lower(barea)){2}.d,3)
    Aover0 = mean(spans0(squeeze(M.(lower(barea)){2}.AB(1,region,i,:))));
    Bover0 = mean(spans0(squeeze(M.(lower(barea)){2}.AB(2,region,i,:))));
    
    biases_time = spans0(squeeze(M.(lower(barea)){2}.dAB(:,region,i,:)));
    dbias_time(i) = nanmean(squeeze(diff(M.(lower(barea)){2}.d([1 5],region,i),[],1)));
    dbias_time_m1(i) = nanmean(squeeze(diff(M.m1{2}.d([1 5],region,i),[],1)));
    AoverB(i) = nanmean(biases_time);
    AoverB_m(i) = nanmean(spans0(squeeze(M.m1{2}.dAB(:,region,i,:))));
    
    AtoB(i) = mean(biases_time(1:round(end/2)))>=1/length(region) & ...
              mean(biases_time(round(end/2):end))<=-1/length(region);
    BtoA(i) = mean(biases_time(1:round(end/2)))<=-1/length(region) & ...
              mean(biases_time(round(end/2):end))>=1/length(region);    
    
    AorB(i) = sum(ismember(biases_time,[1 -1]))/length(region);
    
    favorA(i) = Aover0 > 0.5 & Bover0 < -0.5;
    favorB(i) = Aover0 < -0.5 & Bover0 > 0.5;
    bothAB(i) = Aover0 > 0.5 & Bover0 > 0.5;
    neitherAB(i) = Aover0 < 0.5 & Aover0 > -0.5 & Bover0 < 0.5 & Bover0 > -0.5;
end
anybiasA = find(AoverB > 1/length(region));
anybiasB = find(AoverB < -1/length(region));
notanybias = find(AoverB <= 1/length(region) & AoverB >= -1/length(region));
strongbiasA = find(AoverB > 1/2);
strongbiasB = find(AoverB < -1/2);
weakbiasA = anybiasA(~ismember(anybiasA,strongbiasA));
weakbiasB = anybiasB(~ismember(anybiasB,strongbiasB));
hugebiasA = find(dbias_time > 4); 
hugebiasB = find(dbias_time < -4);
vacillateA = find(AtoB);
vacillateB = find(BtoA);

%% find all biases for all trial blocks
spans0 = @(x) (x(:,1).*x(:,2) > 0).*sign(mean(x,2));
region = mem;
barea = 'pmd';
biases = [];
reps = [];
for j = 1:length(alldays)
    [Aover0,Bover0,AoverB] = deal(zeros(size(M.(lower(barea)){j}.d,3),1));
    for i = 1:size(M.(lower(barea)){j}.d,3)
        Aover0(i) = nanmean(spans0(squeeze(M.(lower(barea)){j}.AB(1,region,i,:))));
        Bover0(i) = nanmean(spans0(squeeze(M.(lower(barea)){j}.AB(2,region,i,:))));

        biases_time = spans0(squeeze(M.(lower(barea)){j}.dAB(:,region,i,:)));
        AoverB(i) = nanmean(biases_time);
    end

    biases.anyA{j}.all = find(AoverB >= 1/length(region));
    biases.anyB{j}.all = find(AoverB <= -1/length(region));
    biases.notany{j}.all = find(AoverB < 1/length(region) & AoverB > -1/length(region));
    biases.strongA{j}.all = find(AoverB > 1/2);
    biases.strongB{j}.all = find(AoverB < -1/2);
    biases.weakA{j}.all = biases.anyA{j}.all(~ismember(biases.anyA{j}.all,biases.strongA{j}.all));
    biases.weakB{j}.all = biases.anyB{j}.all(~ismember(biases.anyB{j}.all,biases.strongB{j}.all));

    reps.anyA{j}.all = find(Aover0 >= 1/length(region));
    reps.anyB{j}.all = find(Bover0 >= 1/length(region));
    reps.AnotB{j}.all = find(Aover0 >= 1/length(region) & Bover0 < 1/length(region));
    reps.BnotA{j}.all = find(Bover0 >= 1/length(region) & Aover0 < 1/length(region));
    reps.AandB{j}.all = find(Aover0 >= 1/length(region) & Bover0 >= 1/length(region));
    reps.AnorB{j}.all = find(Aover0 < 1/length(region) & Bover0 < 1/length(region));
    for k = 1:8
        biases.anyA{j}.targs{k,:} = find(AoverB >= 1/length(region) & targids{j}==k);
        biases.anyB{j}.targs{k,:} = find(AoverB <= -1/length(region) & targids{j}==k);
        biases.notany{j}.targs{k,:} = find(AoverB < 1/length(region) & AoverB > -1/length(region) & targids{j}==k);
        biases.strongA{j}.targs{k,:} = find(AoverB > 1/2 & targids{j}==k);
        biases.strongB{j}.targs{k,:} = find(AoverB < -1/2 & targids{j}==k);
        biases.weakA{j}.targs{k,:} = biases.anyA{j}.targs{k}(~ismember(biases.anyA{j}.targs{k},biases.strongA{j}.targs{k}));
        biases.weakB{j}.targs{k,:} = biases.anyB{j}.targs{k}(~ismember(biases.anyB{j}.targs{k},biases.strongB{j}.targs{k}));

        reps.anyA{j}.targs{k,:} = find(Aover0 >= 1/length(region) & targids{j}==k);
        reps.anyB{j}.targs{k,:} = find(Bover0 >= 1/length(region) & targids{j}==k);
        reps.AnotB{j}.targs{k,:} = find(Aover0 >= 1/length(region) & Bover0 < 1/length(region) & targids{j}==k);
        reps.BnotA{j}.targs{k,:} = find(Bover0 >= 1/length(region) & Aover0 < 1/length(region) & targids{j}==k);
        reps.AandB{j}.targs{k,:} = find(Aover0 >= 1/length(region) & Bover0 >= 1/length(region) & targids{j}==k);
        reps.AnorB{j}.targs{k,:} = find(Aover0 < 1/length(region) & Bover0 < 1/length(region) & targids{j}==k);
    end
    
end
%%
axbias = [];
for ax = 1:4
    axbias.weakA{ax} = [weakbiasA(targids{2}(weakbiasA)==ax); -weakbiasB(targids{2}(weakbiasB)==(ax+4))];
    axbias.strongA{ax} = [strongbiasA(targids{2}(strongbiasA)==ax); -strongbiasB(targids{2}(strongbiasB)==(ax+4))];
    axbias.weakB{ax} = [weakbiasB(targids{2}(weakbiasB)==ax); -weakbiasA(targids{2}(weakbiasA)==(ax+4))];
    axbias.strongB{ax} = [strongbiasB(targids{2}(strongbiasB)==ax); -strongbiasA(targids{2}(strongbiasA)==(ax+4))];
    axbias.notany{ax} = notanybias(ismember(targids{2}(notanybias),[ax ax+4]));
end

%%
reaction_times = [];
postrace = cell(length(alldays),1);
[initang,initerr,finerr] = deal(cell(length(alldays),1));
for i = 1:length(alldays)
    reaction_times.all{i} = alldays(i).tt(:,20)-alldays(i).tt(:,9);
    reaction_times.strongA{i} = reaction_times.all{i}(biases.strongA{i}.all);
    reaction_times.weakA{i} = reaction_times.all{i}(biases.weakA{i}.all);
    reaction_times.strongB{i} = reaction_times.all{i}(biases.strongB{i}.all);
    reaction_times.weakB{i} = reaction_times.all{i}(biases.weakB{i}.all);
    reaction_times.anyA{i} = reaction_times.all{i}(biases.anyA{i}.all);
    reaction_times.anyB{i} = reaction_times.all{i}(biases.anyB{i}.all);
    reaction_times.notany{i} = reaction_times.all{i}(biases.notany{i}.all);
    
    for j = 1:size(alldays(i).tt,1)
        postrace{i}{j,:} = alldays(1).kin.pos(alldays(1).kin.pos(:,1)>alldays(i).tt(j,9) &...
                                alldays(1).kin.pos(:,1)<alldays(i).tt(j,10),2:3);
                            
        inittrace = alldays(1).kin.pos(alldays(1).kin.pos(:,1)>alldays(i).tt(j,9) &...
                                alldays(1).kin.pos(:,1)<alldays(i).tt(j,20),2:3);
        
        if ~isempty(inittrace)
            initang{i}(j,:) = atan2(diff(inittrace([1 end],2),[],1),diff(inittrace([1 end],1),[],1));
        else
            initang{i}(j,:) = NaN;
        end
        
        clc; fprintf('%d - %d/%d\n',i,j,size(alldays(i).tt,1));
    end
        
    finerr{i} = circ_dist(alldays(i).tt(:,17),round(alldays(i).tt(:,17)./(pi/4)).*(pi/4));
    initerr{i} = circ_dist(initang{i},alldays(i).tt(:,17));
    
end
%%
[inittargkap,fintargkap] = deal(zeros(8,1));
for i = 1:8
    inittargkap(i,:) = circ_kappa(initerr{1}(targids{1}==i));
    fintargkap(i,:) = circ_kappa(initerr{1}(targids{1}==i));
end
%%
figure; hold on; 
for i = 1:length(postrace{2})
    if ismember(i,biases.anyA{2})
        plot(postrace{2}{i}(:,1),postrace{2}{i}(:,2),'b'); 
    elseif ismember(i,biases.anyB{2})
        plot(postrace{2}{i}(:,1),postrace{2}{i}(:,2),'r');
    end
end

%%
Comps = [];
for i = 1:8
    opi = mod(i+3,8)+1;
    Comps.pmd{1}(i,:) = nanmean(1./-M.pmd{1}.d(1,:,targids{1}==opi),3);
    Comps.pmd{2}(i,:) = nanmean(1./-M.pmd{1}.d(5,:,targids{1}==i),3);
    
    Comps.m1{1}(i,:) = nanmean(1./-M.m1{1}.d(1,:,targids{1}==opi),3);
    Comps.m1{2}(i,:) = nanmean(1./-M.m1{1}.d(5,:,targids{1}==i),3); 
    
%     Comps.pmd{1} = repmat(min(Comps.pmd{1},[],2),1,size(Comps.pmd{1},2));
%     Comps.pmd{2} = repmat(max(Comps.pmd{2},[],2),1,size(Comps.pmd{2},2));
%     Comps.m1{1} = repmat(min(Comps.m1{1},[],2),1,size(Comps.m1{1},2));
%     Comps.m1{2} = repmat(max(Comps.m1{2},[],2),1,size(Comps.m1{2},2));
    
end

for i = 1:length(M.pmd)
    for j = 1:size(M.pmd{i}.d,3)
        tbtcomp.pmd{i}{1}(:,:,j) = circshift(Comps.pmd{1},5-targids{i}(j),1);
        tbtcomp.pmd{i}{2}(:,:,j) = circshift(Comps.pmd{2},5-targids{i}(j),1);
        
        tbtcomp.m1{i}{1}(:,:,j) = circshift(Comps.m1{1},5-targids{i}(j),1);
        tbtcomp.m1{i}{2}(:,:,j) = circshift(Comps.m1{2},5-targids{i}(j),1);
    end
    
    Mcomp.pmd{i}.d = (1./-M.pmd{i}.d - tbtcomp.pmd{i}{1})./(tbtcomp.pmd{i}{2}-tbtcomp.pmd{i}{1});
    Mcomp.m1{i}.d = (1./-M.m1{i}.d - tbtcomp.m1{i}{1})./(tbtcomp.m1{i}{2}-tbtcomp.m1{i}{1});

    M.pmd{i}.d = 1./-M.pmd{i}.d;
    M.m1{i}.d = 1./-M.m1{i}.d;
%     M.pmd{i}.d = Mcomp.pmd{i}.d;
%     M.m1{i}.d = Mcomp.m1{i}.d;
end

%%
maxspeed = cell(length(alldays),1);
for i = 1:length(alldays)
    for j = 1:size(alldays(i).tt,1)
    
    maxspeed{i}(j,:) = max(sqrt(sum((alldays(1).kin.vel(alldays(1).kin.pos(:,1)>alldays(i).tt(j,9) &...
                                alldays(1).kin.vel(:,1)<alldays(i).tt(j,10),2:3)).^2,2)));
    clc; fprintf('%d: %d/%d\n',i,j,size(alldays(i).tt,1));
    end
end


%%
[var1targ,var2targ,mean1targ,mean2targ,var1targm1,var2targm1,mean1targm1,mean2targm1,mean3targ,mean4targ] = deal(cell(8,1));
for i = 1:8
    var1targ{i} = nanvar(diff(M.pmd{1}.d([1 5],:,targids{1}==i),[],1),[],3); 
    var2targ{i} = nanvar(diff(M.pmd{2}.d([1 5],:,targids{2}==i),[],1),[],3); 
    
    mean1targ{i} = nanmean(diff(M.pmd{1}.d([1 5],:,targids{1}==i),[],1),3);
    mean2targ{i} = nanmean(diff(M.pmd{2}.d([1 5],:,targids{2}==i),[],1),3);
    
    var1targm1{i} = nanvar(diff(M.m1{1}.d([1 5],:,targids{1}==i),[],1),[],3); 
    var2targm1{i} = nanvar(diff(M.m1{2}.d([1 5],:,targids{2}==i),[],1),[],3); 
    
    mean1targm1{i} = nanmean(diff(M.m1{1}.d([1 5],:,targids{1}==i),[],1),3);
    mean2targm1{i} = nanmean(diff(M.m1{2}.d([1 5],:,targids{2}==i),[],1),3);
    
    mean3targ{i} = nanmean(diff(M.pmd{1}.d([1 5],:,targids{3}==i),[],1),3);
    mean4targ{i} = nanmean(diff(M.pmd{1}.d([1 5],:,targids{4}==i),[],1),3);

end  
var1 = cell2mat(var1targ);
var2 = cell2mat(var2targ);
var1m1 = cell2mat(var1targm1);
var2m1 = cell2mat(var2targm1);

mean1 = cell2mat(mean1targ);
mean2 = cell2mat(mean2targ);
mean1m1 = cell2mat(mean1targm1);
mean2m1 = cell2mat(mean2targm1);

maxmean1 = repmat(max(mean1,[],2),1,size(mean1,2));
maxmean1m1 = repmat(max(mean1m1,[],2),1,size(mean1,2));

fano1 = var1./maxmean1;
fano2 = var2./maxmean1;

fano1m1 = var1m1./maxmean1m1;
fano2m1 = var2m1./maxmean1m1;

percchange_pmd = (var2-var1)./var1;
percchange_m1 = (var2m1-var1m1)./var1m1;

pchange_pmd(1,:) = nanmean(percchange_pmd);
[pchange_pmd(2,:),pchange_pmd(3,:)] = boot_bounds(1000,@nanmean,percchange_pmd,2.5,97.5);

pchange_m1(1,:) = nanmean(percchange_m1);
[pchange_m1(2,:),pchange_m1(3,:)] = boot_bounds(1000,@nanmean,percchange_m1,2.5,97.5);

figure; hold on; plot(pchange_pmd(1,:)); plot(pchange_m1(1,:)); 
plot_bounds(pchange_pmd(2,:),pchange_pmd(3,:));
plot_bounds(pchange_m1(2,:),pchange_m1(3,:));



%%
decodedir = cell(length(M.pmd),1);
for blc = 1:length(M.pmd)
    for i = 1:size(M.pmd{blc}.d,3)
        for t = 1:size(M.pmd{blc}.d,2)
            if ~isnan(M.pmd{blc}.d(:,t,i))
                decodedir{blc}(i,t) = find(M.pmd{blc}.d(:,t,i)==max(M.pmd{blc}.d(:,t,i)),1,'first');
            else
                decodedir{blc}(i,t) = nan;
            end
        end
    end
end

%%
[dirdif_1T] = deal(zeros(8,size(M.pmd{1}.d,2)));
dirdif_2T = [];
dfrom1 = [];
for i = 1:8
    dirdif_1T(i,:) = squeeze(nanmean(diff(M.pmd{1}.d([1 5],:,targids{1}==i),[],1),3));
    dirdif_2T.anyA(i,:) = squeeze(nanmean(diff(M.pmd{2}.d([1 5],:,biases.anyA{2}.targs{i}),[],1),3));
    dirdif_2T.anyB(i,:) = squeeze(nanmean(diff(M.pmd{2}.d([1 5],:,biases.anyB{2}.targs{i}),[],1),3));
    dirdif_2T.strongA(i,:) = squeeze(nanmean(diff(M.pmd{2}.d([1 5],:,biases.strongA{2}.targs{i}),[],1),3));
    dirdif_2T.strongB(i,:) = squeeze(nanmean(diff(M.pmd{2}.d([1 5],:,biases.strongB{2}.targs{i}),[],1),3));
    dirdif_2T.weakA(i,:) = squeeze(nanmean(diff(M.pmd{2}.d([1 5],:,biases.weakA{2}.targs{i}),[],1),3));
    dirdif_2T.weakB(i,:) = squeeze(nanmean(diff(M.pmd{2}.d([1 5],:,biases.weakB{2}.targs{i}),[],1),3));
    dirdif_2T.notany(i,:) = squeeze(nanmean(diff(M.pmd{2}.d([1 5],:,biases.notany{2}.targs{i}),[],1),3));
    
    dfrom1.anyA{i} = reshape(squeeze(diff(M.pmd{2}.d([1 5],:,biases.anyA{2}.targs{i}),[],1))',[],size(M.pmd{2}.d,2)) ...
                    -repmat(dirdif_1T(i,:),length(biases.anyA{2}.targs{i}),1);
    dfrom1.anyB{i} = reshape(squeeze(diff(M.pmd{2}.d([1 5],:,biases.anyB{2}.targs{i}),[],1))',[],size(M.pmd{2}.d,2)) ...
                    -repmat(dirdif_1T(i,:),length(biases.anyB{2}.targs{i}),1);
    dfrom1.strongA{i} = reshape(squeeze(diff(M.pmd{2}.d([1 5],:,biases.strongA{2}.targs{i}),[],1))',[],size(M.pmd{2}.d,2))...
                    -repmat(dirdif_1T(i,:),length(biases.strongA{2}.targs{i}),1);
    dfrom1.strongB{i} = reshape(squeeze(diff(M.pmd{2}.d([1 5],:,biases.strongB{2}.targs{i}),[],1))',[],size(M.pmd{2}.d,2))...
                    -repmat(dirdif_1T(i,:),length(biases.strongB{2}.targs{i}),1);
    dfrom1.weakA{i} = reshape(squeeze(diff(M.pmd{2}.d([1 5],:,biases.weakA{2}.targs{i}),[],1))',[],size(M.pmd{2}.d,2))...
                    -repmat(dirdif_1T(i,:),length(biases.weakA{2}.targs{i}),1);
    dfrom1.weakB{i} = reshape(squeeze(diff(M.pmd{2}.d([1 5],:,biases.weakB{2}.targs{i}),[],1))',[],size(M.pmd{2}.d,2))...
                    -repmat(dirdif_1T(i,:),length(biases.weakB{2}.targs{i}),1);
    dfrom1.notany{i} = reshape(squeeze(diff(M.pmd{2}.d([1 5],:,biases.notany{2}.targs{i}),[],1))',[],size(M.pmd{2}.d,2))...
                    -repmat(dirdif_1T(i,:),length(biases.notany{2}.targs{i}),1);
    
end
%%
figure; hold on; 
plot(nanmean(cell2mat(dfrom1.strongA')));
plot(nanmean(cell2mat(dfrom1.strongB')));
plot(nanmean(cell2mat(dfrom1.notany')));

