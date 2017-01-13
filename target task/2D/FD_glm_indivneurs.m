linds = find(AD(2).tt(:,3)==max(AD(2).tt(:,3)));
hinds = find(AD(2).tt(:,3)==min(AD(2).tt(:,3)));

mu_l = circ_mean(AD(2).tt(linds,9));
k_l = circ_kappa(AD(2).tt(linds,9));

sub_h = hinds(rand(length(hinds),1) < circ_vmpdf(AD(2).tt(hinds,9),mu_l,k_l));

random_order_linds = randperm(length(linds));
sub_l = linds(random_order_linds(1:length(sub_h)));

alldays(2).tt = [AD(2).tt(sub_h,:) ; AD(2).tt(sub_l,:)];
Activity_script2;

fd = firing_diffs{1};

% predictors
% distance from prior mean
dfp = abs(circ_dist(alldays(2).tt(:,9),circ_mean(alldays(2).tt(:,9))));

% uncertainty condition
uc = alldays(2).tt(:,3)==min(alldays(2).tt(:,3));

% kinematics (time to target)
ttt = alldays(2).tt(:,7) - alldays(2).tt(:,6);

X1 = [uc];
% % % y = FD(:,20);
% % % 
% % % [b1,dev1,stats1] = glmfit(X,y,'normal');
% % % 
X2 = [dfp];
% % % [b2,dev2,stats2] = glmfit(X2,y,'normal');
% % % 
% % % p_submodel = 1 - chi2cdf(dev2-dev1,1);
X3 = [uc dfp];

%%
[b1,b2,b3,stats1,stats2,stats3] = deal(cell(size(fd,2),length(alldays(1).PMd_units)));
[dev1,dev2,dev3] = deal(zeros(size(fd,2),length(alldays(1).PMd_units)));

[p_un, p_dfp] = deal(nan(size(fd,2),length(alldays(1).PMd_units)));
for j = 1:size(fd,2)
    
    yneurs = fd{j};
    
    for k = 1:size(yneurs,2)
        
        clc; fprintf('time bin: %d/%d\nneuron: %d/%d\n',j,size(fd,2),k,size(yneurs,2));
        y = yneurs(:,k);
    
        [b1{j,k},dev1(j,k),stats1{j,k}] = glmfit(X1,y,'normal');
        [b2{j,k},dev2(j,k),stats2{j,k}] = glmfit(X2,y,'normal');
        [b3{j,k},dev3(j,k),stats3{j,k}] = glmfit(X3,y,'normal');
        
        p_un(j,k) = stats1{j,k}.p(2);
        p_dfp(j,k) = stats2{j,k}.p(2);
        
        b_un(j,k) = b3{j,k}(2);
        b_dfp(j,k) = b3{j,k}(3);
        
        p_submodel(j,k) = 1 - chi2cdf(dev3(j,k)-dev1(j,k),1);
    
    end
%         p_submodel(j) = 1 - chi2cdf(dev2(j)-dev1(j),1);
end
%%
   
p_un_sig = p_un < 0.005;
p_dfp_sig = p_dfp < 0.005;

tot_un_sig = sum(p_un_sig,2);
tot_dfp_sig = sum(p_dfp_sig,2);
    
% const_term = cellfun(@(x) x(1), b1);
% uncert_term = cellfun(@(x) x(2), b1);
% dfp_term = cellfun(@(x) x(3), b1);
% 
% pval_unc = cellfun(@(x) x.p(2), stats1);
% pval_dfp = cellfun(@(x) x.p(3), stats1);

