FIRING = horzcat(firing_absolute{2}{:});%horzcat(Firings{:});
FIRING_CO = horzcat(firing_absolute{1}{:});

% unit_cell = struct2cell(vertcat(neurons{end}{:}));    
% tuning_array = horzcat(unit_cell{1,:})';

unit_cell = struct2cell(vertcat(neurons{end}{:}));    
tuning_array = vertcat(unit_cell{1,:});
% 
modulations = (max(tuning_array,[],2) - min(tuning_array,[],2))';
% 
% minmax_bin = cellfun(@(x) [min(x,[],1);max(x,[],1)],FIRING_CO,'UniformOutput',0);
% minmax_tot = [min(vertcat(minmax_bin{:}),[],1); max(vertcat(minmax_bin{:}),[],1)];

cofire = vertcat(FIRING_CO{:});

goods = find(modulations > MODthresh & mean(cofire)>1);
%goods = find(mean(tuning_array,2)~=inf);

% goods = excites(modulations(excites) > MODthresh);
% goods = inhibs(modulations(inhibs) > MODthresh);

% base_cell = struct2cell(vertcat(neurons{1}{:}));
% baselines = nanmean(horzcat(base_cell{1,:}),1);

%% Likelihoods
likes = flipud(unique(alldays(2).tt(:,3)));
for i = 1:length(likes)
    like_ind{1}{i} = find(alldays(2).tt(:,3)==likes(i));
end
%% Variables

%PDlate = vertcat(PD1{end});
dPD = dfrompd(:,goods);%;{1}{end};

[fire_all,dirs_all,Fire_dir,FIRE_dir] = deal(cell(length(likes),1));

% minmax_bin = cellfun(@(x) [min(x,[],2) max(x,[],2)],TUNES,'UniformOutput',0);
% minmax_tot = [min(horzcat(minmax_bin{:}),[],2), max(horzcat(minmax_bin{:}),[],2)]';

% avrate_tot = nanmean(vertcat(FIRING{:}),1);

spatial_cents = -pi:pi/32:pi;
space_size = diff(spatial_cents(1:2));
[VMparams,VMparamsL,VMparamsH,curveboot,fire_trial,trlCurve,trlBase,trlGain,trlWidth,...
    outL,outH,outT,BM,GM,WM,B,G,W,C] = ...
    deal(cell(length(likes),1));
%% Run
[pfits,curve] = deal(cell(length(FIRING),1));
[base,gain,fwhm] = deal(zeros(size(FIRING{1},1),length(FIRING)));
for bin = 1:length(FIRING)
    
    for i = 1:size(FIRING{1},1)

        FIRE_norm_good = (FIRING{bin}(i,goods))./...
                var(cofire(:,goods));

        [pfits{bin}(i,:),thetas,curve{bin}(i,:),base(i,bin),gain(i,bin),fwhm(i,bin)] = VM_fit_constrained(dPD(i,:),FIRE_norm_good,-pi:0.01:pi);
%         [pfits{bin}(i,:),thetas,curve{bin}(i,:),base(i,bin),gain(i,bin),fwhm(i,bin)] = VM_fit2(dPD(i,:),FIRE_norm_good,-pi:0.01:pi);

    clc; fprintf('fitting...\nbin: %d/%d\ntrial: %d/%d\n',bin,length(FIRING),i,size(FIRING{1},1));
    end
end

%%
figure; hold on; 
cop = {'b','r'};
subplot(1,3,1); hold on; 
for i = 1:length(likes)
    inds = find(alldays(2).tt(:,3)==likes(i));
    bases(i,:) = mean(base(inds,:));
    bases_se(i,:) = std(base(inds,:))./sqrt(size(base,1));
    plot(time_cents,bases(i,:),cop{i});
    patch([time_cents,fliplr(time_cents)],[bases(i,:)+1.96*bases_se(i,:) fliplr(bases(i,:)-1.96*bases_se(i,:))],...
        cop{i},'FaceAlpha',0.25,'EdgeAlpha',0.25);
    
end
subplot(1,3,2); hold on; 
for i = 1:length(likes)
    inds = find(alldays(2).tt(:,3)==likes(i));
    
    gains(i,:) = mean(gain(inds,:));
    gains_se(i,:) = std(gain(inds,:))./sqrt(size(gain,1));
    plot(time_cents,gains(i,:),cop{i});
    patch([time_cents,fliplr(time_cents)],[gains(i,:)+1.96*gains_se(i,:) fliplr(gains(i,:)-1.96*gains_se(i,:))],...
        cop{i},'FaceAlpha',0.25,'EdgeAlpha',0.25);
    
end
subplot(1,3,3); hold on; 
for i = 1:length(likes)
    inds = find(alldays(2).tt(:,3)==likes(i));
    fwhms(i,:) = mean(fwhm(inds,:));
    fwhms_se(i,:) = std(fwhm(inds,:))./sqrt(size(fwhm,1));
    plot(time_cents,fwhms(i,:),cop{i});
    patch([time_cents,fliplr(time_cents)],[fwhms(i,:)+1.96*fwhms_se(i,:) fliplr(fwhms(i,:)-1.96*fwhms_se(i,:))],...
        cop{i},'FaceAlpha',0.25,'EdgeAlpha',0.25);
    
    cvs{i} = cellfun(@(x) mean(x(inds,:)),curve,'UniformOutput',0);
    
    avcurves{i} = vertcat(cvs{i}{:});
    
end
minc = min([min(min(avcurves{1})) min(min(avcurves{2}))]);
maxc = max([max(max(avcurves{1})) max(max(avcurves{2}))]);

figure; hold on; 
subplot(1,2,1); 
imagesc(time_cents,spatial_cents,avcurves{1}',[minc,maxc]); 
subplot(1,2,2); 
imagesc(time_cents,spatial_cents,avcurves{2}',[minc,maxc]);

