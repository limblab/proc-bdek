% %
aligns = {'target','go'};
timebins = {[-800:50:800], [0:50:1000]};
arrays = {'PMd','M1'};

cutility = cell(length(arrays),length(aligns));

for barea = 1:length(arrays)
for epoch = 1:length(aligns)

brain_area = arrays{barea};
time_bins = timebins{epoch};
time_align = aligns{epoch};
prediction_day_indices = [2];
 
unit_counts = cell(length(prediction_day_indices),1);
unit_rasts = cell(length(prediction_day_indices),1);
for z = 1:length(prediction_day_indices)

    prediction_day = prediction_day_indices(z);
    day_units = eval(sprintf('alldays(%d).%s_units',prediction_day,brain_area));
    if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; M1_ids(end-1:end) = []; end

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
%margins = zeros(size(tbt{z},2),length(tbt));
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
            
            prior_prob = length(inds)/(num_trials-1);
            
            final_prob{binnum}(tr,i) = prod(10*probs_given_spikes)*prior_prob;
            neuron_prob{binnum}(tr,i,:) = 10*probs_given_spikes*prior_prob;
          
        end
        
        Induct_class = Like_uncs(final_prob{binnum}(tr,:)==max(final_prob{binnum}(tr,:)));
        
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
        
        temp_marg = sqrt(sum(neuron_trans_prob{binnum}(tr,:,:).^2,2));
        discrim = dot(neuron_trans_prob{binnum}(tr,:,:),ones(size(neuron_trans_prob{binnum}(tr,:,:))));
        
        cor_discrim = dot(neuron_trans_prob{binnum}(tr,:,:),repmat(correct_ind,[1 1 size(neuron_trans_prob{binnum},3)]));
        
        correct_discrim{binnum}(tr,:) = ...
            dot(neuron_trans_prob{binnum}(tr,:,:)./repmat(temp_marg,[1 size(neuron_trans_prob{binnum},2) 1]),...
                repmat(correct_ind,[1 1 size(neuron_trans_prob{binnum},3)]));
        
        margins_i{binnum}(tr,:) = temp_marg(:)';
        un_discrim{binnum}(tr,:) = acos(discrim(:)'./(temp_marg(:)' * norm(ones(1,length(Like_uncs)))));
        %correct_discrim2{binnum}(tr,:) = pi/2 - acos(cor_discrim(:)'./(temp_marg(:)'));
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
        
        pref_of_i{binnum}(:,i) = final_prob{binnum}(:,i) == max(final_prob{binnum},[],2);
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
fprintf('Done\n');
%%
cur_tt = alldays(prediction_day_indices(z)).tt;
Like_uncs = unique(cur_tt(:,3));
num_trials = size(cur_tt,1);
% Cycle through trials

final_probk = zeros(size(tbt{z},2),length(Like_uncs));
final_probstd = zeros(size(tbt{z},2),length(Like_uncs));
pref_of_ki = zeros(size(tbt{z},2),length(Like_uncs));
pref_of_std = zeros(size(tbt{z},2),length(Like_uncs));
actual_ki = zeros(size(tbt{z},2),length(Like_uncs));
un_discrimstd = zeros(size(tbt{z},2),1);
for tr = 1:size(tbt{z},2)

    clc;
    fprintf('bin: %d/%d\ntrial: %d/%d\n',binnum,length(tbt),tr,size(tbt{z},2));

    % X endpoint train
    Xk = cur_tt(:,10);
    % X time to target
    Xttt = cur_tt(:,7)-cur_tt(:,6);
    % X slice std
    Xstd = cur_tt(:,11);
    
    % X test
    Xktest = cur_tt(tr,10);
    Xttttest = cur_tt(tr,7) - cur_tt(tr,6);
    Xstdtest = cur_tt(tr,11);
    
    for i = 1:length(Like_uncs)
        inds = find(cur_tt(:,3)==Like_uncs(i));
        inds(inds==tr) = [];
        
        training_data_m = mean(Xk(inds));
        training_data_k = circ_kappa(Xk(inds));
        
        training_ttt_m = mean(Xttt(inds));
        training_ttt_std = std(Xttt(inds));
        
        training_std_m = mean(Xstd(inds));
        training_std_std = std(Xstd(inds));
        
        probs_given_kin = circ_vmpdf(Xktest,training_data_m,training_data_k);
        probs_given_ttt = normpdf(Xttttest,training_ttt_m,training_ttt_std);
        probs_given_std = normpdf(Xstdtest,training_std_m,training_std_std);
        
        prior_probk = length(inds)/num_trials;

        final_probk(tr,i) = probs_given_kin * probs_given_ttt * prior_probk;
        final_probstd(tr,i) = probs_given_std * prior_probk;
    end
    
    temp_margstd = norm(final_probstd(tr,:));
    discrimstd = dot(final_probstd(tr,:),ones(1,length(Like_uncs)));
    un_discrimstd(tr) = acos(discrimstd/(temp_margstd * norm(ones(1,length(Like_uncs)))));
end
    
for i = 1:length(Like_uncs)

    pref_of_ki(:,i) = final_probk(:,i) == max(final_probk,[],2);
    pref_of_std(:,i) = final_probstd(:,i) == max(final_probstd,[],2);
    actual_ki(:,i) = cur_tt(:,3) == Like_uncs(i);
    
end

comparison = pref_of_ki + actual_ki;
comparisonstd = pref_of_std + actual_ki;
correct_trialsk = find(max(comparison,[],2)==2);
correct_trialsstd = find(max(comparisonstd,[],2)==2);

perc_correctk = length(correct_trialsk)./size(cur_tt,1);
perc_correctstd = length(correct_trialsstd)./size(cur_tt,1);

%%

neur_util = zeros(length(unit_rasts{1}),length(tbt));
sess_util = zeros(length(tbt),size(tbt{1},2));
neur_cutil = zeros(length(unit_rasts{1}),length(tbt));
sess_cutil = zeros(length(tbt),size(tbt{1},2));
for i = 1:length(tbt)
    neur_util(:,i) = mean(un_discrim{i})';
    neur_cutil(:,i) = mean(correct_discrim{i})';
end
for i = 1:size(tbt{1},2)
    for j = 1:length(tbt)
        sess_util(j,i) = mean(un_discrim{j}(i,:));
        sess_cutil(j,i) = mean(correct_discrim{j}(i,:));
    end
end

cutility{barea,epoch} = neur_cutil;
end
end

maxsize = max(max(cell2mat(cutility)));
minsize = dot([1 1 1]/norm([1 1 1]),[0 0 1]);
%minsize = min(min(cell2mat(cutility)));
rangesize = maxsize - minsize;

orig_place = cd;
denloc = cell(length(arrays),length(aligns));

size_maps = cell(length(arrays),length(aligns));
for barea = 1:length(arrays)
    for epoch = 1:length(aligns)

        brain_area = arrays{barea};
        time_bins = timebins{epoch};
        time_align = aligns{epoch};
        
        mkdir(sprintf('C:\\Users\\limblab\\Desktop\\figures\\%s\\%s\\movie\\%s\\%s\\',monkey,modaye,brain_area,time_align));
        for tt = 1:size(cutility{barea,epoch},2)
            
            map_input = 40*(cutility{barea,epoch}(:,tt) - minsize)/rangesize + 1;
            
            [size_maps{barea,epoch}{tt}, ~, ~, ~] = ArrayMapPlot('Mihili',brain_area,eval(sprintf('%s_ids',brain_area)),...
                map_input,1,ones(length(eval(sprintf('%s_ids',brain_area))),1),'',1); 
            cd((sprintf('C:\\Users\\limblab\\Desktop\\figures\\%s\\%s\\movie\\%s\\%s\\',monkey,modaye,brain_area,time_align)));
            saveas(fig_hand,sprintf('%s_%s_%s_%d.png',monkey,brain_area,time_align,tt),'png');
            close(fig_hand);
            cd(orig_place);
            clear fig_hand;
            clc; fprintf('%d/%d\n',tt,size(cutility{barea,epoch},2));

%             thismap = size_maps{barea,epoch}{tt}; thismap(thismap < 0) = 0;
%             summap = sum(thismap,3);
%             xsum = sum(summap)/sum(sum(summap));
%             ysum = sum(summap,2)/sum(sum(summap));
%             
%             xhighdens = mean(xsum.*[1:10]);
%             yhighdens = mean(ysum.*[10:-1:1]');
%             
%             denloc{barea,epoch}(tt,:) = [xhighdens, yhighdens];
        end
    end
end
cd(orig_place);