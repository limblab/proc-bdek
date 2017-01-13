FIRING = horzcat(firing_absolute{2}{:});%horzcat(Firings{:});
FIRING_CO = horzcat(firing_absolute{1}{:});

unit_cell = struct2cell(vertcat(neurons{end}{:}));    
%tuning_array = vertcat(unit_cell{1,:});
tuning_array = horzcat(unit_cell{1,:})';

modulations = (max(tuning_array,[],2) - min(tuning_array,[],2))';

cofire = vertcat(FIRING_CO{:});

goods = find(modulations > MODthresh & mean(cofire)>1);

%% Likelihoods
likes = flipud(unique(alldays(2).tt(:,3)));
for i = 1:length(likes)
    like_ind{1}{i} = find(alldays(2).tt(:,3)==likes(i));
end
%% Variables

dPD = dfrompd(:,goods);

[fire_all,dirs_all,Fire_dir,FIRE_dir,fire_trial] = deal(cell(length(likes),1));

spatial_cents = -pi:pi/16:pi;
space_size = diff(spatial_cents(1:2));
%% Run
for lik = 1:length(likes)
    
    inds = find(alldays(2).tt(:,3)==likes(lik));
    dirs_all{lik} = reshape(dPD(inds,:),[],1);
    for bin = 1:length(FIRING)

        FIRE_norm_good = (FIRING{bin}(inds,goods))./...
                    repmat(var(cofire(:,goods)),length(inds),1);

        fire_all{lik}(:,bin) = reshape(FIRE_norm_good,[],1);
        
        fire_trial{lik}{bin} = FIRE_norm_good;

    end
    
    for dir = 1:length(spatial_cents)
        
        dirinds = find(abs(circ_dist(dirs_all{lik},spatial_cents(dir)))<(space_size/2));
        Fire_dir{lik}(dir,:) = nanmean(fire_all{lik}(dirinds,:),1);
        FIRE_dir{lik}{dir} = fire_all{lik}(dirinds,:);
    end 
end

%%
[outE,outL,B,G,W] = deal(cell(length(likes),1));
[PARAMS_E,PARAMS_L] = deal(nan(3,length(likes)));
for lik = 1:length(likes)
    good_dirs = find(~isnan(sum(Fire_dir{lik},2)));
    
%     xsfireE = repmat(spatial_cents(good_dirs)',4,1);
%     xsfireL = repmat(spatial_cents(good_dirs)',4,1);
%     ysfireE = reshape(Fire_dir{lik}((good_dirs),1:4),[],1);
%     ysfireL = reshape(Fire_dir{lik}((good_dirs),6:9),[],1);
    
    xsfireE = spatial_cents(good_dirs)';
    xsfireL = spatial_cents(good_dirs)';
    ysfireE = Fire_dir{lik}(good_dirs,1);
    ysfireL = Fire_dir{lik}(good_dirs,3);
    
    [~,~,~,PARAMS_E(1,lik),PARAMS_E(2,lik),PARAMS_E(3,lik)] = VM_fit3(xsfireE,ysfireE,-pi:0.01:pi);
    [outE{lik}(1,:),outE{lik}(2,:)] = boot_bounds(100,@VM_helper_func_kap,[xsfireE, ysfireE],2.5,97.5);

    [~,~,~,PARAMS_L(1,lik),PARAMS_L(2,lik),PARAMS_L(3,lik)] = VM_fit3(xsfireL,ysfireL,-pi:0.01:pi);
    [outL{lik}(1,:),outL{lik}(2,:)] = boot_bounds(100,@VM_helper_func_kap,[xsfireL, ysfireL],2.5,97.5);     
    
end

for lik = 1:length(likes)
    
    B{lik} = [PARAMS_E(1,lik) PARAMS_L(1,lik); outE{lik}(:,1) outL{lik}(:,1)];
    G{lik} = [PARAMS_E(2,lik) PARAMS_L(2,lik); outE{lik}(:,2) outL{lik}(:,2)];
    W{lik} = [PARAMS_E(3,lik) PARAMS_L(3,lik); outE{lik}(:,3) outL{lik}(:,3)];
    
end
%%
[PARAMS] = cell(length(likes),1);
for lik = 1:length(likes)
    
    for i = 1:length(FIRING)
        ysfire = Fire_dir{lik}(:,i);

        [~,~,~,PARAMS{lik}(1,i),PARAMS{lik}(2,i),PARAMS{lik}(3,i)] = VM_fit3(spatial_cents',ysfire,-pi:0.01:pi);
    end
end