% figure; hold on; 

% for i = 1:length(minds)
%     
%     mi = minds(i);
%     plot(6*cos(0:0.01:2*pi),6*sin(0:0.01:2*pi),'k');
%     plot(7*cos(0:0.01:2*pi),7*sin(0:0.01:2*pi),'k');
%     for j = 1:5
%         plot([6 7]*cos(alldays(2).slices(mi,j)),[6 7]*sin(alldays(2).slices(mi,j)),'b','LineWidth',3);
%     end
%     plot([6 7]*cos(alldays(2).tt(mi,9)),[6 7]*sin(alldays(2).tt(mi,9)),'r','LineWidth',3);
%     
%     plot(6.5*cos(alldays(2).tt(mi,10)),6.5*sin(alldays(2).tt(mi,10)),'g.','MarkerSize',15);
% 
%     pause; cla
% end

%%
% figure; hold on; 

thets = 0.01:0.01:2*pi;
kern_K = 2.^(0:9);
convest = zeros(length(alldays(2).tt),length(kern_K));
avcon = zeros(length(alldays(2).tt),1);
convd = cell(length(alldays(2).tt),1);
avcrves = zeros(length(alldays(2).tt),length(thets));
for i = 1:length(alldays(2).tt)
    
    slcs = mod(alldays(2).slices(i,:),2*pi);
    if isnan(sum(slcs))
        convest(i,:) = nan(1,length(kern_K));
        convd{i} = nan(length(kern_K),length(thets));
    else
    
        dels = zeros(length(thets),1); dels(round(100*slcs)) = 1;

        for k = 1:length(kern_K)
            kern = circ_vmpdf(thets,0,kern_K(k));

            convd{i}(k,:) = cconv(dels,kern,length(kern));
            maxloc = thets(convd{i}(k,:)==max(convd{i}(k,:)));
            convest(i,k) = maxloc(1);
            
            np = findpeaks(convd{i}(k,:));
            if np == 5
                convd{i}(k,:) = NaN;
            end
        end
        avcrve = mean(convd{i});
        avcrves(i,:) = avcrve./sum(avcrve).*length(avcrve);
        maxavcrve = thets(avcrve==max(avcrve));
        avcon(i,:) = maxavcrve(1);
        
%         plot(dels); 
%         plot(convd','b');
%         plot(avcrve,'k','LineWidth',3);
%         plot([100 100]*alldays(2).tt(i,10),[0 max(convd(:))],'k','LineWidth',5);
%         pause; cla
    end
end

%%
errf = @(x,y) nansum(circ_dist(x,y).^2)./sum(~isnan(x));
like_inds = {linds, minds, hinds};

contender_list = [convest avcon alldays(2).tt(:,9)];
for i = 1:length(kern_K)
    
    for lik = 1:length(like_inds)
        
        averr(i,lik) = errf(contender_list(like_inds{lik},i),alldays(2).tt(like_inds{lik},10));
    end
end

%%
figure; hold on; 
for i = 1:size(avcrves,1)
    minma = [min(avcrves(i,:)) max(avcrves(i,:))];
    plot(thets,avcrves(i,:)); 
    plot([1 1]*alldays(2).tt(i,10),minma,'r','LineWidth',3); 
    plot([1 1]*avcon(i),minma,'b','LineWidth',3);
    pause; 
    cla; 
end

%%
trials2use = [linds;minds;hinds];


ws = ones(size(convd{1},1),1);

weightsum = @(w,X) thets(find(nansum(repmat(w,1,size(X,2)).*X)==max(nansum(repmat(w,1,size(X,2)).*X)),1,'first'));

wsfun = @(w) nansum(circ_dist(cellfun(@(x) weightsum(w,x),convd(trials2use)),alldays(2).tt(trials2use,10)).^2)./sum(~isnan(trials2use));

[Ws,minval] = fminsearch(wsfun,ws);

%%
patrd = zeros(length(alldays(2).tt),1);
top20 = zeros(length(alldays(2).tt),125);
for i = 1:length(alldays(2).tt)
    
    rdi = round(100*alldays(2).tt(i,10));
    
    patrd(i,:) = avcrves(i,rdi);
    
    [~,likelis] = sortrows(avcrves(i,:)');
    
    top20(i,:) = thets(likelis(end-124:end));
    
end
%%
[pk1est,pkclosest,pknum,pkperc,pkrange] = deal(zeros(length(alldays(2).tt),1));
pks = cell(length(alldays(2).tt),1);
for i = 1:length(alldays(2).tt)
    [pkval,peaklocs] = findpeaks(avcrves(i,:));
    val_loc = flipud(sortrows([pkval' peaklocs']));
   
    
    maxpk = find(abs(circ_dist(thets(val_loc(:,2)),alldays(2).tt(i,10)))==...
                    min(abs(circ_dist(thets(val_loc(:,2)),alldays(2).tt(i,10)))));
    
    if ~isempty(peaklocs)
        pknum(i) = maxpk;
        pkclosest(i) = thets(val_loc(pknum(i),2));
        pk1est(i) = thets(val_loc(1,2));
        pkperc(i) = val_loc(pknum(i),1)./max(val_loc(:,1));
        pks{i} = val_loc;
        
        pkcombs = combnk(deg2rad(pks{i}(:,2)),2);
        if ~isempty(pkcombs)
            pkrange(i) = max(abs(circ_dist(pkcombs(:,1),pkcombs(:,2))));
        else
            pkrange(i) = 0;
        end
    else
        pkclosest(i) = NaN;
        pkperc(i) = NaN;
        pk1est(i) = NaN;
        pknum(i) = NaN;
        pks{i} = nan;
        pkrange(i) = nan;
    end

end
    
%% Get motor var with dir
dirs = unique(alldays(1).tt(:,2));
[kap_motor, errsqr] = deal(zeros(length(dirs),1));
for i = 1:length(dirs)
    dirinds = find(alldays(1).tt(:,2)==dirs(i));
    kap_motor(i) = circ_kappa(alldays(1).tt(dirinds,10));
    errsqr(i) = circ_mean(circ_dist(alldays(1).tt(dirinds,10),alldays(1).tt(dirinds,2)).^2);
end


%% fit kernel convolution combined with prior distrib. 
prior_dist = circ_vmpdf(thets,circ_mean(alldays(2).tt(:,2)),circ_kappa(alldays(2).tt(:,2))); 
last_trial = [NaN; alldays(2).tt(1:end-1,2)];
prev_error = [NaN; circ_dist(alldays(2).tt(1:end-1,10),alldays(2).tt(1:end-1,2))];

conv_prior = cell(size(avcrves,1),1);
for i = 1:size(avcrves,1)
    conv_prior{i} = [avcrves(i,:); prior_dist'];
end

trials2use = hinds;%1:size(alldays(2).tt,1);
pmu = circ_mean(alldays(2).tt(:,2));
pk = circ_kappa(alldays(2).tt(:,2));


ws = ones(size(conv_prior{1},1),1);

weightsum = @(w,X) thets(find(nansum(repmat(w,1,size(X,2)).*X)==max(nansum(repmat(w,1,size(X,2)).*X)),1,'first'));
postfunc = @(mu,k,lik) thets(find(circ_vmpdf(thets,mu,k).*lik'==max(circ_vmpdf(thets,mu,k).*lik'),1,'first'));

wsfun = @(w) nansum(circ_dist(cellfun(@(x) weightsum(w,x),conv_prior(trials2use)),alldays(2).tt(trials2use,10)).^2)./sum(~isnan(trials2use));

kfun = @(k) nansum(circ_dist(cellfun(@(x) postfunc(pmu,k,x(1,:)),conv_prior(trials2use)),alldays(2).tt(trials2use,10)).^2)./sum(~isnan(trials2use));


[Ws,minval] = fminsearch(wsfun,ws);
[Ks,minval_k] = fminsearch(kfun,pk);

pred_loc_sum = cellfun(@(X) weightsum(Ws,X),conv_prior);
pred_loc_prod = cellfun(@(X) postfunc(pmu,Ks,X(1,:)),conv_prior);

%% Choose option closest to prior
for i = 1:size(avcrves)
    
    peak2prior = abs(circ_dist(pks{i}(:,2)./100,circ_mean(alldays(2).tt(:,2))));
    
    best_peak(i) = pks{i}(peak2prior==min(peak2prior),2)./100;
end



