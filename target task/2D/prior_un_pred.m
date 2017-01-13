%%
brain_area = 'PMd';
time_bins = [0 800];
time_align = 'target';
prediction_day_indices = [2 3];
 
unit_counts = cell(length(prediction_day_indices),1);
unit_rasts = cell(length(prediction_day_indices),1);
for z = 1:length(prediction_day_indices)

    prediction_day = prediction_day_indices(z);
    day_units = eval(sprintf('alldays(%d).%s_units',prediction_day,brain_area));
    if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end

    day_pred = prediction_day;
    tt_pred = alldays(day_pred).tt;

    % find rasters for each unit
    unit_pred = cell(length(day_units),1);
    for q = 1:length(day_units)
        
        clc;
        fprintf('Day: %d/%d\nUnit: %d/%d\n',z,length(prediction_day_indices),q,length(day_units));

        [rast_out,rast_inds] = raster_plot(day_units{q},tt_pred,[time_bins(1)./1000 time_bins(end)./1000],...
            time_align,0,'none');

        out = nan(size(vertcat(rast_out{:})));
        for m = 1:length(rast_out)
            out(rast_inds{m},:) = rast_out{m};
        end
        out(isnan(out))=0;
        rast = zeros(size(out,1),length(time_bins)-1);
        binsizes = time_bins - time_bins(1); binsizes(1) = 1;
        for v = 1:length(time_bins)-1
            rast(:,v) = sum(out(:,binsizes(v):binsizes(v+1)),2);
        end

        unit_rasts{z}{q} = rast;
        %unit_counts{z}{q} = rast./(diff(time_bins)/1000);
    end
end
%% Population Prediction of Prior condition
all_trials = vertcat(unit_rasts{:});
comb_rasts = cell(length(unit_rasts));
for i = 1:length(all_trials)
    comb_rasts{i} = vertcat(all_trials{:,i});
end

tbt = cell(1,length(time_bins)-1);
for i = 1:size(comb_rasts{1},2)
    for j = 1:length(comb_rasts)
        tbt{i}(j,:) = comb_rasts{j}(:,i)';
    end
end

prinums = [];
hacky = [1 2];
for i = 1:length(unit_rasts)
    prinums = [prinums; hacky(i)*ones(size(unit_rasts{i}{1},1),1)];
end
% Cycle through trials
final_prob = cell(1,length(tbt));
kinem_prob = cell(1,length(tbt));
neuron_prob = cell(1,length(tbt));
num_trials = length(prinums);
margins = zeros(size(tbt{1},2),length(tbt));
full_tt = vertcat(alldays(prediction_day_indices).tt);

for binnum = 1:length(tbt)
    for tr = 1:size(tbt{1},2)
        
        clc;
        fprintf('bin: %d/%d\ntrial: %d/%d\n',binnum,length(tbt),tr,size(tbt{1},2));

        % X train
        X = tbt{binnum};
        Xk = full_tt(:,10);
        Xttt = full_tt(:,7)-full_tt(:,6);
        
        % X test
        Xtest = tbt{binnum}(:,tr);
        Xk_test = full_tt(tr,10);
        Xttt_test = full_tt(tr,7)-full_tt(tr,6);
        
        for i = 1:length(unit_rasts)
            inds = find(prinums==hacky(i));
            inds(inds==tr) = [];
            
            training_data = mean(X(:,inds),2);
            
            training_k_m = mean(Xk(inds));
            training_k_k = circ_kappa(Xk(inds));

            training_ttt_m = mean(Xttt(inds));
            training_ttt_std = std(Xttt(inds));

            probs_given_kin = circ_vmpdf(Xk_test,training_k_m,training_k_k);
            probs_given_ttt = normpdf(Xttt_test,training_ttt_m,training_ttt_std);
            
            probs_given_spikes = poisspdf(Xtest,training_data);
            prior_prob = length(inds)/num_trials;
            
            final_prob{binnum}(tr,i) = prod(probs_given_spikes)*prior_prob;
            kinem_prob{binnum}(tr,i) = prod(probs_given_kin)*prod(probs_given_ttt)*prior_prob;
        end
    end
end
    
pref_of_i = cell(1,length(tbt));
pref_kin = cell(1,length(tbt));
actual_i = cell(1,length(tbt));
perc_correct = zeros(1,length(tbt));
prefi = cell(1,length(tbt));
acti = cell(1,length(tbt));
p_corr = zeros(size(X,1),length(tbt));
correct_trials = cell(length(tbt),1);
correct_trials_k = cell(length(tbt),1);
for binnum = 1:length(tbt)
    
    for i = 1:length(unit_rasts)
        
        pref_of_i{binnum}(:,i) = final_prob{binnum}(:,i) == max(final_prob{binnum},[],2);
        pref_kin{binnum}(:,i) = kinem_prob{binnum}(:,i) == max(kinem_prob{binnum},[],2);
       
        actual_i{binnum}(:,i) = prinums == i;
       
    end
    
    comparison = pref_of_i{binnum} + actual_i{binnum};
    correct_trials{binnum} = find(max(comparison,[],2)==2);
    
    perc_func = @(x) sum(x)/(length(x));
    
    perc_correct(binnum) = perc_func(correct_trials{binnum});
    
    comparison_k = pref_kin{binnum} + actual_i{binnum};
    correct_trials_k{binnum} = find(max(comparison_k,[],2)==2);
    
    perc_correctk(binnum) = perc_func(correct_trials_k{binnum});
 
    [pc_L(binnum), pc_H(binnum)] = boot_bounds(10000,perc_func,correct_trials{binnum},2.5,97.5);
    
end
   
%%
numbins = 10;
within_bin = length(prinums)/numbins;
n = zeros(length(tbt),numbins);
b = n;
nk = n; 
bk = n;
for i = 1:length(tbt)
    [n(i,:),b(i,:)] = hist(correct_trials{i},numbins);
    [nk(i,:),bk(i,:)]=hist(correct_trials_k{i},numbins);
end
figure; hold on; 
plot(b',n'./within_bin,'.-');
xtimebins = .5*(time_bins(1:end-1)+time_bins(2:end));
legend('Hold','Delay');
hold on; 

switch_prior = find(diff(prinums)~=0);
plot(repmat(switch_prior,1,2)',repmat([0 1.1],length(switch_prior),1)','k--');
ylim([0 1.1]);

xlabel('Trials','FontSize',16);
ylabel('% Correct','FontSize',16);
title(sprintf('%s',brain_area),'FontSize',14);

figure; hold on; 
plot(bk',100*((n-nk)./nk)','.-');
legend('Hold','Delay');
hold on; 
plot(repmat(switch_prior,1,2)',repmat([-100 100],length(switch_prior),1)','k--');
ylim([-100 120]);
xlabel('Trials','FontSize',16);
ylabel('% Improvement','FontSize',16);
title(sprintf('Improvement over kinematic (%s)',brain_area),'FontSize',14);

figure; hold on;
plot(bk',nk'./within_bin,'.-');
%%
% maxind = zeros(size(p_corr,1),1);
% maxpred = zeros(size(p_corr,1),1);
% maxdelay = zeros(size(p_corr,1),1);
% for i= 1:size(p_corr,1)
%     
%     maxpred(i) = max(p_corr(i,:));
%     maxind(i) = mean(xtimebins(p_corr(i,:)==maxpred(i)));
%     maxdelay(i) = max(p_corr(i,xtimebins>0));
% end
% 
% hold_period_perf = reshape(p_corr(:,xtimebins > 0),sum(xtimebins > 0)*size(p_corr,1),1);
% hp_mean = mean(hold_period_perf);
% hp_se = 1.96*std(hold_period_perf);
% 
% hp_thresh = hp_mean + hp_se;
% over_thresh = sum(maxdelay > hp_thresh)/size(p_corr,1); 
% 
% perc_units_better = sum(maxdelay > perc_correctk)/length(day_units);
