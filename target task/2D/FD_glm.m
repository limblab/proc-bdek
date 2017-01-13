
pdi = 2;
FD = zeros(length(firing_diffs{1}{1}),length(firing_diffs{1}));
for i = 1:length(firing_diffs{1})
    FD(:,i) = nanmean(firing_diffs{1}{i},2);
end

% Output variable (deviation from expected firing
FD_all = reshape(FD,[],1);

% predictors
% distance from prior mean
dfp = abs(circ_dist(alldays(pdi).tt(:,9),circ_mean(alldays(pdi).tt(:,9))));
dfp_rep = repmat(dfp,size(FD,2),1);

% uncertainty condition
uc = alldays(pdi).tt(:,3)==min(alldays(pdi).tt(:,3));
uc_rep = repmat(uc,size(FD,2),1);

% kinematics (time to target)
ttt = alldays(pdi).tt(:,7) - alldays(pdi).tt(:,6);
ttt_rep = repmat(ttt,size(FD,2),1);

% Trial number
tn = (1:size(alldays(pdi).tt,1))';
tn_rep = repmat(tn,size(FD,2),1);

% time within trial (SEPARATE VARIABLES)
twt = zeros(length(FD_all),size(FD,2));
for i = size(FD,2)
    twt(:,i) = zeros(length(FD_all),1);
    twt((size(FD,1)*(i-1)+1):(size(FD,1)*i),i) = 1;
end
    
% 
% X = [dfp_rep uc_rep ttt_rep tn_rep];
% y = FD_all;

X1 = [uc dfp];
% % % y = FD(:,20);
% % % 
% % % [b1,dev1,stats1] = glmfit(X,y,'normal');
% % % 
X2 = [uc];
% % % [b2,dev2,stats2] = glmfit(X2,y,'normal');
% % % 
% % % p_submodel = 1 - chi2cdf(dev2-dev1,1);
X3 = [dfp];

XC = [uc dfp uc.*dfp];

[b1,b2,b3,bc,bb,bbint,bbstats,bbsub,bbsubint,rsub,bbsubstats,r,stats1,...
    stats2,stats3,statsc,bcr,b1r,b2r,bbsub2,bbsubint2,rsub2,bbsubstats2] = deal(cell(size(FD,2),1));
[dev1,dev2,dev3,devc,devpart2] = deal(zeros(size(FD,2),1));
LI = like_ind{1};

[devfull,devpart,fitimprove_unc,fitimprove_dfp] = deal(zeros(size(FD,2),1000));
[fitimproveLH_unc, fitimproveLH_dfp, IMPROVE_un,IMPROVE_dfp] = deal(zeros(size(FD,2),2));
for j = 1:size(FD,2)
    clc; fprintf('%d/%d\n',j,size(FD,2));
    y = FD(:,j);
    
    [b1{j},dev1(j),stats1{j}] = glmfit(X1,y,'normal');
    [b2{j},dev2(j),stats2{j}] = glmfit(X2,y,'normal');
    [b3{j},dev3(j),stats3{j}] = glmfit(X3,y,'normal');
    [bc{j},devc(j),statsc{j}] = glmfit(XC,y,'normal');
    [bb{j},bbint{j},r{j},~,bbstats{j}] = regress(y,[ones(size(XC,1),1) XC]);
    [bbsub{j},bbsubint{j},rsub{j},~,bbsubstats{j}] = regress(y,[ones(size(X2,1),1) X2]);
    [bbsub2{j},bbsubint2{j},rsub2{j},~,bbsubstats2{j}] = regress(y,[ones(size(X3,1),1) X3]);
    
%     [~,samp] = bootstrp(1000,[],1:length(uc));
%     for brep = 1:size(samp,2)
%         [bcr{j}(:,brep), devfull(j,brep)] = glmfit(XC(samp(:,brep),:),y(samp(:,brep)),'normal');
%         [b1r{j}(:,brep), devpart(j,brep)] = glmfit(X3(samp(:,brep),:),y(samp(:,brep)),'normal'); %dfp alone
%         [b2r{j}(:,brep), devpart2(j,brep)]= glmfit(X2(samp(:,brep),:),y(samp(:,brep)),'normal'); %unc alone
%         
%         fitimprove_unc(j,brep) = (devpart(j,brep) - devfull(j,brep))./devfull(j,brep);
%         fitimprove_dfp(j,brep) = (devpart2(j,brep) - devfull(j,brep))./devfull(j,brep);
%     end
    

    
    R2imp_un = rsub2{j}.^2 - r{j}.^2;
    R2imp_dfp = rsub{j}.^2 - r{j}.^2;
    
    [IMPROVE_un(j,1), IMPROVE_un(j,2)] = boot_bounds(1000,@sum,R2imp_un,2.5,97.5);
    [IMPROVE_dfp(j,1), IMPROVE_dfp(j,2)] = boot_bounds(1000,@sum,R2imp_dfp,2.5,97.5);
        
%     sorted_unc = sortrows(fitimprove_unc(j,:)'); 
%     sorted_dfp = sortrows(fitimprove_dfp(j,:)');
%     
%     fitimproveLH_unc(j,1) = sorted_unc(round(0.025*1000)); fitimproveLH_unc(j,2) = sorted_unc(round(0.975*1000));
%     fitimproveLH_dfp(j,1) = sorted_dfp(round(0.025*1000)); fitimproveLH_dfp(j,2) = sorted_dfp(round(0.975*1000));
%     
%     p_submodel(j) = 1 - chi2cdf(dev2(j)-dev1(j),1);

end
% %%
% for i = 1:length(LI)
%     
%     likeinds = LI{i};
% 
%     for j = 1:size(FD,2)
%         y = FD(likeinds,j);
%         Xin = X3(likeinds,:);
%         [blike{i,j},devlike{i,j},statslike{i,j}] = glmfit(Xin,y,'normal');
%     end
%     
%     B_const(i,:) = cellfun(@(x) x(1),blike(i,:));
%     B_dfp(i,:) = cellfun(@(x) x(2),blike(i,:));
%     
%     B_dfpSE(i,:) = cellfun(@(x) x.se(2),statslike(i,:));
% 
% end
% 
% 
%     %%
% const_term = cellfun(@(x) x(1), b1);
% uncert_term = cellfun(@(x) x(2), b1);
% dfp_term = cellfun(@(x) x(3), b1);
% 
% pval_unc = cellfun(@(x) x.p(2), stats1);
% pval_dfp = cellfun(@(x) x.p(3), stats1);

%%


