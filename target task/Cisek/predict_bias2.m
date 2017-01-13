%-%
brain_area = 'PMd';
time_bins = [500 1000];
time_align = 6;
Train_block = 7;
Test_block = 9;
 
%-% 
day_units = eval(sprintf('alldays(1).%s_units',brain_area));
if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end

tt_train = alldays(Train_block).tt;
tt_test = alldays(Test_block).tt;
% Find training rasters for each unit, separated by target location 
unit_pred = cell(length(day_units),1);
train_rast = cell(8,length(day_units));
test_rast = cell(1,length(day_units));
for q = 1:length(day_units)
    clc; fprintf('Unit: %d/%d\n',q,length(day_units));
    [rast_out,rast_inds] = raster_plotCK(day_units{q},tt_train,[time_bins(1)./1000 time_bins(end)./1000],...
                           time_align,0,'none',[],13); % Group into correct target locations
    
    [rast_outT,rast_indsT] = raster_plotCK(day_units{q},tt_test,[time_bins(1)./1000 time_bins(end)./1000],...
                           time_align,0,'none',[],13); % Group into correct target locations       
    ro = vertcat(rast_outT{:});
    ri = vertcat(rast_indsT{:});
    [~,ranked] = sortrows(ri);
    
    test_rast{q} = ro(ranked,:);                   
    train_rast(:,q) = rast_out;
end
targ_locs = unique(tt_train(:,13));
find_targ_ind = @(x) find(abs(circ_dist(x,targ_locs))==min(abs(circ_dist(x,targ_locs))));
%%
%-% Loop through test trials and categorize direction
[P1,P2,Pref_prob] = deal(zeros(size(tt_test,1),length(time_bins)-1));
[class,classrat,classind,codeoff,targrat,scalingfact,scaled_rat] = deal(zeros(size(tt_test,1),length(time_bins)-1)); 
indivrat = nan(size(tt_test,1),length(alldays(1).PMd_units));
[valds,IRATS,dfts] = deal(zeros(size(tt_test,1),length(day_units)));
[IP1,IP2] = deal(cell(size(tt_test,1),1));
for i = 1:size(tt_test,1)
    clc; fprintf('Trial: %d/%d (%.0f%%)\n ',i,size(tt_test,1),100*i/size(tt_test,1));
    p1 = find_targ_ind(tt_test(i,17)); p2 = find_targ_ind(tt_test(i,17)+pi);
    %   p1 = find_targ_ind(tt_test(i,2)); p2 = find_targ_ind(tt_test(i,3));
    p1S = mod([p1-2 p1],8)+1; p2S = mod([p2-2 p2],8)+1;
    
    Train1 = train_rast(p1,:); Train2 = train_rast(p2,:);
    Train1S = train_rast(p1S,:); Train2S = train_rast(p2S,:);
    
    [probs1,probs2] = deal(zeros(length(day_units),length(time_bins)-1));
    
    training = train_rast([p1,p2],:);
    %training = train_rast([p1,p1S,p2,p2S],:);
    training_bin = cellfun(@(x) bin_array(x,size(x,1),length(time_bins)-1,'sum'),training,'UniformOutput',0);
    %training_bin = cellfun(@(x) bin_array(x,size(x,1),1,'sum'),training,'UniformOutput',0);
    
    train_means = cell2mat(cellfun(@(x) nanmean(x,1)',training_bin,'UniformOutput',0));
    train_sum = cell2mat(cellfun(@(x) nansum(x,1)',training_bin,'UniformOutput',0));
    trialnumbs = repmat([size(Train1{1},1);size(Train2{1},1)],1,length(day_units));
    %train_means = train_sum./((length(time_bins)-1)*trialnumbs);
    
    valid_neurs = find(prod(train_means,1)~=0);
    
    dftaxis = abs(dot([cos(PD1{1}) sin(PD1{1})],...
        repmat([cos(tt_test(i,17)) sin(tt_test(i,17))],length(day_units),1),2));
    
    dftaxis_rep = repmat(dftaxis',size(train_means,1)/2,1);
    
    dftaxis_val = dftaxis(valid_neurs);
    
    dfts(i,:) = dftaxis';
%     
    tmeans1 = train_means(1:end/2,:);
    tmeans2 = train_means((end/2+1):end,:);
%     tmeans1 = repmat(train_means(1,:),length(time_bins)-1,1);
%     tmeans2 = repmat(train_means(2,:),length(time_bins)-1,1);
    
    testing_bin = cellfun(@(x) bin_array(x(i,:),1,length(time_bins)-1,'sum')',test_rast,'UniformOutput',0);
    test_in = cell2mat(testing_bin);
    
    indiv_probs1 = poisspdf(test_in,tmeans1);
    indiv_probs2 = poisspdf(test_in,tmeans2);
    
    IP1{i} = indiv_probs1;
    IP2{i} = indiv_probs2;
    
    probs_given_train1 = prod(5.*indiv_probs1(:,valid_neurs),2);
    probs_given_train2 = prod(5.*indiv_probs2(:,valid_neurs),2);
    
    i_rats = indiv_probs1(:,valid_neurs)./indiv_probs2(:,valid_neurs);
    i_rats_dnorm = i_rats;%.^(dftaxis_rep(:,valid_neurs));
    
    IRATS(i,:) = log(indiv_probs1./indiv_probs2);
    valds(i,valid_neurs) = 1;
    %helper_func = @(k) -log(prod(i_rats.^(k*dftaxis_rep(:,valid_neurs))));
    %kap(i) = fmincon(helper_func,1,[],[],[],[],0.1);
        
    tot_rat = prod(i_rats_dnorm,2);
    lograt = log(tot_rat);
    
    class(i,:) = ((lograt < 0) + 1)';
    %classrat(i,:) = ((c2./c1).^(-2.*class(i,:)' + 3))';
    targrat(i,:) = lograt';
    
    codeoff(i,:) = class(i,:)*pi - pi;
    for j = 1:length(time_bins)-1
        classind(i,j) = find_targ_ind(tt_test(i,17)+codeoff(i,j));
    end
    %fprintf(sprintf('%.1f/%.1f/%.1f\n',a,b,c));
%     for q = 1:length(day_units)
%         
%             test_vec = bin_array(test_rast{q}(i,:),1,length(time_bins)-1,'sum');
%             train_vec = [bin_array(train_rast{p1,q},1,length(time_bins)-1,'sum')/size(train_rast{p1,q},1);...
%                          bin_array(train_rast{p2,q},1,length(time_bins)-1,'sum')/size(train_rast{p1,q},1)];
%                      
%             probs1(q,:) = poisspdf(test_vec,train_vec(1,:));
%             probs2(q,:) = poisspdf(test_vec,train_vec(2,:));         
%     end
%     P1(i,:) = nanmean(probs1);
%     P2(i,:) = nanmean(probs2);
%     Pref_prob(i,:) = P1(i,:)-P2(i,:);
end
%%
% tP1 = mean(P1,2);
% tP2 = mean(P2,2);
% 
% pref1 = tP1./tP2;
% 
% pref = pi-pi*(pref1>1);
% figure; hold on; 
% colrs = {'b','r'};
% for i = 1:size(indivrat,1)
%     [~,yp] = linehist(-1:.1:1,indivrat(i,:),colrs{cue1or2(i)}); 
%     plot(medneur(i)*[1 1],[0 max([yp yp2])],'k','LineWidth',5);
%     pause; cla; 
% end
%     


%helper_func = @(k) -mean(sum(log(IRATSg).*(dftsg>k),2)>0);


