% %%
brain_area = 'all';
time_bins = [-200:100:600];
time_align = 'target';
prediction_day_indices = [2];
 
unit_counts = cell(length(prediction_day_indices),1);
unit_rasts = cell(length(prediction_day_indices),1);
for z = 1:length(prediction_day_indices)

    prediction_day = prediction_day_indices(z);
    day_units = eval(sprintf('alldays(1).%s_units',brain_area));
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
%% Population Prediction of Likelihood condition

z = 1;
tbt = cell(1,length(time_bins)-1);
for i = 1:size(unit_rasts{z}{1},2)
    for j = 1:length(unit_rasts{z})
        tbt{i}(j,:) = unit_rasts{z}{j}(:,i)';
    end
end

cur_tt = alldays(prediction_day_indices(z)).tt;
Like_uncs = unique(cur_tt(:,3));
num_trials = size(cur_tt,1);
% Cycle through trials
final_prob = cell(1,length(tbt));
neuron_prob = cell(1,length(tbt));
trans_prob = cell(1,length(tbt));
neuron_trans_prob = cell(1,length(tbt));
margins_i = cell(1,length(tbt));
un_discrim = cell(1,length(tbt));
correct_discrim = cell(1,length(tbt));
class_margin = cell(1,length(tbt));
neur_margin = cell(1,length(tbt));
allprev_prob = cell(1,length(tbt));
prev_pri = cell(1,length(tbt));
neur_correct = cell(1,length(tbt));
%margins = zeros(size(tbt{z},2),length(tbt));
nprob = cell(1,length(tbt));

Like_list = alldays(prediction_day_indices).tt(:,3);
for i = 1:length(Like_uncs)
    likeconds{i} = find(Like_list == Like_uncs(i));
end

for binnum = 1:length(tbt)
    for tr = 1:size(tbt{z},2)
        
        clc;
        fprintf('bin: %d/%d\ntrial: %d/%d\n',binnum,length(tbt),tr,size(tbt{z},2));

        % X train
        X = tbt{binnum};
        bad_uns = find(sum(X,2)==0);
        
        % X test
        Xtest = tbt{binnum}(:,tr);
        
        for i = 1:length(Like_uncs)
            inds = find(cur_tt(:,3)==Like_uncs(i));
            inds(inds==tr) = [];
            training_data = mean(X(:,inds),2);
            probs_given_spikes = poisspdf(Xtest,training_data);
            probs_given_spikes(probs_given_spikes==0) = 1;
            
            prior_prob = 0.5; %length(inds)/(num_trials-1);
            
            final_prob{binnum}(tr,i) = prod(10*probs_given_spikes)*prior_prob;
            neuron_prob{binnum}(tr,i,:) = 10*probs_given_spikes*prior_prob;
            nprob{binnum}(tr,i,:) = 10*probs_given_spikes;
          
        end
        
        nonemfp = final_prob(1:binnum);
        all_probs = cellfun(@(x) x(tr,:),nonemfp,'UniformOutput',0);
        prev_probs = vertcat(all_probs{1:binnum});
        prev_probs = prev_probs./repmat(max(prev_probs,[],2),1,length(Like_uncs));
        prev_pri{binnum}(tr,:) = prod(prev_probs,1);
        
        allprev_prob{binnum}(tr,:) = prev_pri{binnum}(tr,:).*final_prob{binnum}(tr,:);
        
        Induct_class = Like_uncs(final_prob{binnum}(tr,:)==max(final_prob{binnum}(tr,:)));
        %Induct_class = Like_uncs(allprev_prob{binnum}(tr,:)==max(allprev_prob{binnum}(tr,:)));
        class_margin{binnum}(tr) = max(final_prob{binnum}(tr,:))./(sum(final_prob{binnum}(tr,:)));
        neur_margin{binnum}(:,tr) = max(neuron_prob{binnum}(tr,:,:),[],2)./sum(neuron_prob{binnum}(tr,:,:),2);

        for i = 1:length(Like_uncs)
            Like_list = cur_tt(:,3);
            Like_list(tr) = Induct_class;
            inds = find(Like_list==Like_uncs(i));
            training_data = mean(X(:,inds),2);
            probs_given_spikes = poisspdf(Xtest,training_data);
            probs_given_spikes(probs_given_spikes==0) = 1;
            prior_prob = length(inds)/num_trials;
            
            trans_prob{binnum}(tr,i) = prod(probs_given_spikes)*prior_prob;
            neuron_trans_prob{binnum}(tr,i,:) = probs_given_spikes*prior_prob;
            
        end
       
        correct_ind = (Like_uncs == cur_tt(tr,3))';
        aa = reshape(neuron_prob{binnum}(tr,correct_ind,:),[],1);
        bb = reshape(sum(neuron_prob{binnum}(tr,:,:),2),[],1);
        neur_correct{binnum}(:,tr) = aa./bb;
        %neur_correct{binnum}(:,tr) = neuron_prob{binnum}(tr,correct_ind,:)./sum(neuron_prob{binnum}(tr,:,:),2);
        
        temp_marg = sqrt(sum(neuron_trans_prob{binnum}(tr,:,:).^2,2));
        discrim = dot(neuron_trans_prob{binnum}(tr,:,:),ones(size(neuron_trans_prob{binnum}(tr,:,:))));
        
        cor_discrim = dot(neuron_trans_prob{binnum}(tr,:,:),repmat(correct_ind,[1 1 size(neuron_trans_prob{binnum},3)]));
        
        margins_i{binnum}(tr,:) = temp_marg(:)';
        un_discrim{binnum}(tr,:) = acos(discrim(:)'./(temp_marg(:)' * norm(ones(1,length(Like_uncs)))));
        correct_discrim{binnum}(tr,:) = pi/2 - acos(cor_discrim(:)'./(temp_marg(:)'));
    end
    
    %margins(:,binnum) = (max(final_prob{binnum},[],2)-min(final_prob{binnum},[],2))./max(final_prob{binnum},[],2);
    %dpp(:,binnum) = 1 - dot(final_prob{binnum},trans_prob{binnum},2);    
    
end
    
pref_of_i = cell(1,length(tbt));
actual_i = cell(1,length(tbt));
perc_correct = zeros(1,length(tbt));
prefi = cell(1,length(tbt));
acti = cell(1,length(tbt));
p_corr = zeros(size(X,1),length(tbt));
boot_perc = zeros(2,length(tbt));
fprintf('Finishing...\n');
for binnum = 1:length(tbt)
    
    for i = 1:length(Like_uncs)
        
        pref_of_i{binnum}(:,i) = allprev_prob{binnum}(:,i) == max(allprev_prob{binnum},[],2);
        %pref_of_i{binnum}(:,i) = final_prob{binnum}(:,i) == max(final_prob{binnum},[],2);
        actual_i{binnum}(:,i) = cur_tt(:,3) == Like_uncs(i);
        
        for j = 1:size(X,1)
            
            prefi{binnum}{j}(:,i) = neuron_prob{binnum}(:,i,j) == max(neuron_prob{binnum}(:,:,j),[],2);
            acti{binnum}{j}(:,i) = cur_tt(:,3) == Like_uncs(i);
            
        end
    end
    
    comparison = pref_of_i{binnum} + actual_i{binnum};
    correct_trials = find(max(comparison,[],2)==2);
    correct_binned = max(comparison,[],2)==2;
    
    perc_func = @(x) sum(x)/length(x);
    
    perc_correct(binnum) = perc_func(correct_binned);
    [boot_perc(1,binnum), boot_perc(2,binnum)] = boot_bounds(10000,perc_func,correct_binned,2.5,97.5);
    
    
    for j = 1:size(X,1)
        
        comp = prefi{binnum}{j} + acti{binnum}{j};
        corrs = find(max(comp,[],2)==2);
        
        p_corr(j,binnum) = length(corrs)./size(cur_tt,1);
    end
    
end
%%
final_neur_prob = nprob{end};
p = prod(final_neur_prob,3);

ratio_n = zeros(size(final_neur_prob,3),1);
for i = 1:size(final_neur_prob,3)
    fnpi = reshape(final_neur_prob(:,:,i),size(final_neur_prob,1),[]);
    
    correct_neurp = sum(fnpi.*acti{end}{end},2);
    incorrect_neurp = sum(fnpi.*(acti{end}{end}==0),2);
    
    av_corr_p = mean(correct_neurp);
    av_incor_p = mean(incorrect_neurp);
    
    %ratio_n(i) = av_corr_p./av_incor_p;
    ratio_n(i) = mean(correct_neurp./incorrect_neurp);

end
%%
perccordn = zeros(size(final_neur_prob,3)-1,1);
removed = [];
for ns = 1:(size(final_neur_prob,3)-1)
    clear pc_wi;
    clc; fprintf('dropped: %d\n',ns);
    for i = 1:size(final_neur_prob,3) % Through all neurons
        
        if ns ==1
            removed = 1;
        else
            removed = [i bestneur(1:(ns-1))];
        end

        % probabilities with dropped neurons
        missingi = final_neur_prob;
        missingi(:,:,removed) = NaN;
        pwithouti = prod(missingi(:,:,~ismember(1:size(missingi,3),removed)),3);

        % classification without included neurons
        choice_woi = pwithouti==repmat(max(pwithouti,[],2),1,length(Like_uncs));
        
        % correct classifications
        correcttrials = max(choice_woi + acti{end}{end},[],2)==2;
        corr_margins(i) = mean(max(pwithouti,[],2) - min(pwithouti,[],2));
        
        for lclass = 1:length(Like_uncs)
            corrclass(lclass) = mean(correcttrials(likeconds{lclass}));
        end
        
        correct_probs = sum(pwithouti.*acti{end}{end},2);
        incorrect_probs = sum(pwithouti.*(acti{end}{end}==0),2);
        
        pratio = mean(correct_probs)./mean(incorrect_probs);

        %pc_wi(i) = mean(correcttrials); % percent correct
        pc_wi(i) = pratio;
    end

    pc_wi(removed) = inf;
    
    worsti = find(pc_wi == min(pc_wi));
    % average margin of correct classification
    if length(worsti) > 1  
        bestmarg = find(corr_margins(worsti) == max(corr_margins(worsti),[],2),1,'last'); 
        bestneur(ns) = worsti(bestmarg);
    else
        bestneur(ns) = worsti;
    end
    
    perccordn(ns) = min(pc_wi);

end
    
%%

fprintf('Done\n');

figure; hold on; 
xtimebins = .5*(time_bins(1:end-1)+time_bins(2:end));
plot(xtimebins,perc_correct,'b','LineWidth',2);
patch([xtimebins fliplr(xtimebins)],[boot_perc(1,:) fliplr(boot_perc(2,:))],'b','FaceAlpha',0.5,'EdgeAlpha',0.5);
xlabel('Time from Target On (ms)','FontSize',16);
ylabel('% Correct Classification','FontSize',16);
title(sprintf('%s',brain_area),'FontSize',14);
ylim([0 1]);

